#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <algorithm>
#include <cmath>

// ===== 你项目已有的头文件 =====
#include "types.hpp"
#include "gf128.hpp"
#include "hash_prg.hpp"
#include "poly.hpp"
#include "rbokvs.hpp"
#include "layout.hpp"
#include "wire.hpp"
#include "lagrange.hpp"

// —— 分阶段计时结构 —— //
struct Timings {
  double s12_ms{}, s13_ms{}, s14_ms{}, s31_ms{}, total_ms{};
};

// *** COMM *** 新增通信统计结构
struct Comm {
    uint64_t S13 = 0;   // OKVS 上传通信量
    uint64_t S14 = 0;   // σ 查询通信量
    uint64_t total() const { return S13 + S14; }
};

// 简便构造 Block128
static inline Block128 X(u64 a, u64 b){ return Block128{a,b}; }

// —— 单次跑完整流程 —— //
static double run_once(
    int n, int S, int k,
    double eps_okvs, uint32_t w, double eps_hash,
    u64 salt_tag,
    Timings* t = nullptr,
    Comm* comm = nullptr   // *** COMM ***
){
  using clk = std::chrono::high_resolution_clock;
  auto g0 = clk::now();   // 开始

  // ====== 构造每方集合 ======
  std::vector<std::vector<Block128>> Xs(n+1);
  int common = std::max(1, S/2);
  std::vector<Block128> common_elems;
  common_elems.reserve(common);
  for(int tt=0; tt<common; ++tt) common_elems.push_back(X(0, 100+tt));

  for(int i=1; i<=n; ++i){
    Xs[i] = common_elems;
    for(int u=0; u<S-common; ++u){
      Xs[i].push_back(X(1000+i, 100000 + 1000*i + u));
    }
  }

  // ====== S12 ======
  std::vector<std::vector<KV>> kv_all(n+1);
  std::vector<std::vector<Tag128>> tag_all(n+1);

  for(int i=1; i<=n; ++i){
    Block128 alpha_i = hash_to_field_u64((u64)i);
    kv_all[i].reserve(Xs[i].size());
    tag_all[i].reserve(Xs[i].size());

    for(const auto& x : Xs[i]){
      constexpr u64 poly_salt = 0xC0FFEEULL;   // 固定全局盐(只用于多项式系数)
      PRG prg(mix_seed(x, poly_salt));         // 只依赖 x
      Block128 fxi = poly_eval_fx_at_i(x, alpha_i, k, prg);
      kv_all[i].push_back({x, fxi});
      tag_all[i].push_back(tag_of(x, salt_tag));
    }
  }
  auto g1 = clk::now();

  // ====== S13: 各方编码 OKVS ======
  std::vector<RBOKVS> okvs(n+1);
  std::vector<size_t> ni(n+1);
  for(int i=1; i<=n; ++i){
    ni[i] = Xs[i].size();

    OKVSParams p;
    p.m = static_cast<size_t>(std::ceil((1.0 + eps_okvs) * ni[i]));
    if (p.m < (size_t)w + 1) p.m = (size_t)w + 1;
    p.w = w;
    p.seed_r1 = 0xA1B2C3D400000000ULL ^ (uint64_t)i;
    p.seed_r2 = 0x0F1E2D3C00000000ULL ^ ((uint64_t)i << 8);

    okvs[i] = RBOKVS::Encode(kv_all[i], p);

    // *** COMM *** S13：每方上传 OKVS 表大小
    if (comm) {
        uint64_t okvs_bytes = okvs[i].byte_size();
      comm->S13 += okvs_bytes;
    }
  }
  auto g2 = clk::now();

  // ====== S14 ======
  size_t M = 0; 
  for(int i=1;i<=n;++i) M = std::max(M, ni[i]);
  size_t B = size_t(eps_hash * M) + 1;

  std::vector<HashTableTi> Ts(n+1, HashTableTi{std::vector<Bucket>(B)});

  for(int i=1;i<=n;++i){
    for(size_t j=0;j<ni[i];++j){
      const Block128& x   = kv_all[i][j].key;
      const Block128& fxi = kv_all[i][j].val;

      std::vector<std::pair<int, Block128>> sigmas;
      sigmas.reserve(n-1);

      for(int g=1; g<=n; ++g){
        if(g==i) continue;

        // *** COMM *** i → g 发送 x（16 字节）
        if (comm) comm->S14 += sizeof(Block128);

        Block128 sig = okvs[g].Decode(x);

        // *** COMM *** g → i 发送 σ（16 字节）
        if (comm) comm->S14 += sizeof(Block128);

        sigmas.emplace_back(g, sig);
      }

      insert_element_Ti(Ts[i], B, n, i, x, tag_all[i][j], fxi, sigmas, salt_tag);
    }
  }
  auto g3 = clk::now();

  // ====== S31–S33（无通信，不统计） ======
  const int agg = n; (void)agg;
  std::set<u64> result_hash;
  auto tag_key = [](const Tag128& t)->u64 { return t.hi ^ (t.lo<<1); };

  for(size_t eta=0; eta<B; ++eta){
    std::vector<Share> pool; pool.reserve(8);
    for(int i=1;i<=n;++i)
      for(const auto &sh : Ts[i].table[eta].items) pool.push_back(sh);

    std::unordered_map<u64, std::vector<Share>> by_tag;
    for(const auto &sh: pool) by_tag[tag_key(sh.tag)].push_back(sh);

    for(auto& kv : by_tag){
      auto &vec = kv.second;

      std::vector<std::pair<Block128,Block128>> pts;
      std::unordered_set<int> seen_party;
      pts.reserve(vec.size());

      for(const auto &sh : vec){
        if(seen_party.insert(sh.party_id).second){
          pts.push_back({ gf128::from_u64((u64)sh.party_id), sh.fx_i });
        }
      }
      if((int)pts.size() < k) continue;

      std::vector<std::pair<Block128,Block128>> use_pts(pts.begin(), pts.begin()+k);

      Block128 s_rec = gf128::lagrange_at_gf128(use_pts, gf128::zero());

      bool ok = true;
      for(const auto &pr : use_pts){
        Block128 check = gf128::lagrange_at_gf128(use_pts, pr.first);
        if ((check.hi ^ pr.second.hi) | (check.lo ^ pr.second.lo)){ ok = false; break; }
      }

      if(ok) result_hash.insert(kv.first);
    }
  }
  auto g4 = clk::now();

  // ====== 写阶段耗时 ======
  double s12 = std::chrono::duration<double,std::milli>(g1-g0).count();
  double s13 = std::chrono::duration<double,std::milli>(g2-g1).count();
  double s14 = std::chrono::duration<double,std::milli>(g3-g2).count();
  double s31 = std::chrono::duration<double,std::milli>(g4-g3).count();
  double total = std::chrono::duration<double,std::milli>(g4-g0).count();
  if (t) { t->s12_ms=s12; t->s13_ms=s13; t->s14_ms=s14; t->s31_ms=s31; t->total_ms=total; }

  return total;
}


int main(){
    // ====== 协议固定参数 ======
    const double eps_okvs = 0.05;
    const uint32_t w = 192;
    const double eps_hash = 1.3;
    const u64 salt_tag = 0x13572468abcdef99ULL;

    // ====== 实验变量 ======
    std::vector<int> m_values = {4, 16, 32, 64, 256, 1024,4096};
    std::vector<int> n_values = {5, 10, 20, 30, 40, 50};
    std::vector<int> t_values = {2, 3, 5, 7};
    const int reps = 5;

    // 百分位函数
    auto percentile = [](std::vector<double> v, double p){
        size_t N = v.size();
        size_t idx = (size_t)std::round((N - 1) * p);
        std::nth_element(v.begin(), v.begin()+idx, v.end());
        return v[idx];
    };

    // ====== 输出运行时间 ======
    std::ofstream out("rt_m_t_n.csv", std::ios::out | std::ios::trunc);
    out << "m_fixed,t_eff,n,rt_median_s,rt_p2p5_s,rt_p97p5_s\n";

    // *** COMM *** 输出通信量
    std::ofstream out_comm("comm_m_t_n.csv", std::ios::out | std::ios::trunc);
    out_comm << "m_fixed,t_eff,n,S13_bytes,S14_bytes,total_bytes\n";

    // ====== 主循环 ======
    for(int m : m_values){
        for(int t : t_values){
            for(int n : n_values){

                std::vector<double> v;
                std::vector<uint64_t> S13s, S14s, Totals;
                v.reserve(reps);
                S13s.reserve(reps);
                S14s.reserve(reps);
                Totals.reserve(reps);

                for(int r=0; r<reps; ++r){
                    Comm comm;   // *** COMM ***
                    double ms =
                        run_once(n, m, t, eps_okvs, w, eps_hash, salt_tag,
                                 nullptr, &comm);

                    v.push_back(ms / 1000.0);
                    S13s.push_back(comm.S13);
                    S14s.push_back(comm.S14);
                    Totals.push_back(comm.total());
                }

                double med  = percentile(v, 0.5);
                double p2   = percentile(v, 0.025);
                double p97  = percentile(v, 0.975);

                // 写入运行时间
                out << m << "," << t << "," << n << ","
                    << med << "," << p2 << "," << p97 << "\n";

                // *** COMM *** 通信量中位数
                auto percentile_u64 = [&](std::vector<uint64_t> v, double p){
                    size_t N = v.size();
                    size_t idx = (size_t)std::round((N - 1) * p);
                    std::nth_element(v.begin(), v.begin()+idx, v.end());
                    return v[idx];
                };

                uint64_t S13_med = percentile_u64(S13s, 0.5);
                uint64_t S14_med = percentile_u64(S14s, 0.5);
                uint64_t Tot_med = percentile_u64(Totals, 0.5);

                out_comm << m << "," << t << "," << n << ","
                         << S13_med << "," << S14_med << "," << Tot_med
                         << "\n";

                std::cout << "[m="<<m<<", t="<<t<<", n="<<n
                          << "] COMM S13="<<S13_med
                          << " S14="<<S14_med
                          << " TOTAL="<<Tot_med << "\n";
                std::cout << "[m="<<m<<", t="<<t<<", n="<<n
          << "] RT median="<<med<<" s"
          << " p2.5="<<p2<<" s"
          << " p97.5="<<p97<<" s"
          << "\n";

            }
        }
    }

    std::cout << "✅ Wrote rt_m_t_n.csv and comm_m_t_n.csv\n";
    return 0;
}

