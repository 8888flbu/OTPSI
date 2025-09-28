#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
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

// 简便构造 Block128
static inline Block128 X(u64 a, u64 b){ return Block128{a,b}; }

// —— 单次跑完整流程：给定 n(参与方), S(每方集合大小)，返回总耗时(ms)
static double run_once(
    int n, int S, int k,
    double eps_okvs, uint32_t w, double eps_hash,
    u64 salt_tag
){
  using clk = std::chrono::high_resolution_clock;
  auto t0 = clk::now();

  // ====== 构造每方集合：S 个元素，约 50% 交集 + 50% 各自独有（可稳定复现）======
  std::vector<std::vector<Block128>> Xs(n+1);
  int common = std::max(1, S/2);
  std::vector<Block128> common_elems;
  common_elems.reserve(common);
  for(int t=0; t<common; ++t) common_elems.push_back(X(0, 100+t)); // 所有方共享

  for(int i=1; i<=n; ++i){
    Xs[i] = common_elems;
    for(int u=0; u<S-common; ++u){
      // 为每个参与方补充各自的独有元素（不会与他人重叠）
      Xs[i].push_back(X(1000+i, 100000 + 1000*i + u));
    }
  }

  // ====== S12: 生成 (x, f_x(i)) 与 tag ======
  std::vector<std::vector<KV>>       kv_all(n+1);
  std::vector<std::vector<Tag128>>   tag_all(n+1);

  for(int i=1; i<=n; ++i){
    Block128 alpha_i = hash_to_field_u64((u64)i);
    kv_all[i].reserve(Xs[i].size());
    tag_all[i].reserve(Xs[i].size());

    for(const auto& x : Xs[i]){
      PRG prg(mix_seed(x, 0xC0FFEEULL ^ (u64)i));
      Block128 fxi = poly_eval_fx_at_i(x, alpha_i, k, prg);
      kv_all[i].push_back({x, fxi});
      tag_all[i].push_back(tag_of(x, salt_tag));
    }
  }

  // ====== S13: 各方编码 OKVS ======
  std::vector<RBOKVS> okvs(n+1);
  std::vector<size_t> ni(n+1);
  for(int i=1; i<=n; ++i){
    ni[i] = Xs[i].size();
    OKVSParams p;
    p.m = static_cast<size_t>(std::ceil((1.0 + eps_okvs) * ni[i]));
    if (p.m < (size_t)w + 1) p.m = (size_t)w + 1;    // 保证 m - w >= 1
    p.w = w;
    p.seed_r1 = 0xA1B2C3D400000000ULL ^ (uint64_t)i;
    p.seed_r2 = 0x0F1E2D3C00000000ULL ^ ((uint64_t)i << 8);

    okvs[i] = RBOKVS::Encode(kv_all[i], p);
  }

  // ====== S14: 解码 σ_{i,j}^γ 并构建各方哈希表 T_i ======
  size_t M = 0; for(int i=1;i<=n;++i) M = std::max(M, ni[i]);
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
        Block128 sig = okvs[g].Decode(x);
        sigmas.emplace_back(g, sig);
      }
      insert_element_Ti(Ts[i], B, n, i, x, tag_all[i][j], fxi, sigmas, salt_tag);
    }
  }

  // ====== S31–S33: 终端聚合方做 GF(2^128) 拉氏插值与一致性验证 ======
  const int agg = n;
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

      // 取前 k 个点做插值（也可改成随机抽 k 个做多次鲁棒恢复）
      std::vector<std::pair<Block128,Block128>> use_pts(pts.begin(), pts.begin()+k);
      Block128 s_rec = gf128::lagrange_at_gf128(use_pts, gf128::zero());

      // 一致性检查
      bool ok = true;
      for(const auto &pr : use_pts){
        Block128 check = gf128::lagrange_at_gf128(use_pts, pr.first);
        if ( (check.hi ^ pr.second.hi) | (check.lo ^ pr.second.lo) ){
          ok = false; break;
        }
      }
      if(ok) result_hash.insert(kv.first);
    }
  }

  auto t1 = clk::now();
  return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

int main(){
  // ====== 固定协议参数（与原代码一致，可按需调整）======
  const int k = 2;                        // 阈值
  const double   eps_okvs = 0.05;         // RB-OKVS 负载参数
  const uint32_t w        = 192;          // RB-OKVS 带宽
  const double   eps_hash = 1.3;          // 哈希表负载因子
  const u64      salt_tag = 0x13572468abcdef99ULL;

  // ====== 横轴（参与方个数）与多条折线的集合大小 ======
  std::vector<int> party_counts = {5, 10, 20, 30, 40, 50};
  std::vector<int> set_sizes    = {4, 16, 32};   // 多条线：4/16/32 elements
  const int reps = 3; // 每个点跑多次取最优（或平均），降低抖动

  std::ofstream out("scaling_bench.csv", std::ios::out | std::ios::trunc);
  if(!out){ std::cerr << "cannot open scaling_bench.csv\n"; return 1; }
  out << "set_size,n_parties,runtime_s\n";

  for(int S : set_sizes){
    for(int n : party_counts){
      double best_ms = 1e100;
      for(int r=0; r<reps; ++r){
        double ms = run_once(n, S, k, eps_okvs, w, eps_hash, salt_tag);
        if(ms < best_ms) best_ms = ms;
      }
      out << S << "," << n << "," << (best_ms/1000.0) << "\n";
      std::cout << "[S="<<S<<", n="<<n<<"] runtime = " << (best_ms/1000.0) << " s\n";
    }
  }
  std::cout << "Wrote scaling_bench.csv\n";
  return 0;
}

