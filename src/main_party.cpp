#include <iostream>
#include <unordered_map>
#include <map>
#include <set>
#include "types.hpp"
#include "gf128.hpp"
#include "hash_prg.hpp"
#include "poly.hpp"
#include "rbokvs.hpp"
#include "layout.hpp"
#include "wire.hpp"

// 简易辅助：把原始“明文集合”构造为几个元素
static Block128 X(u64 a, u64 b){ return Block128{a,b}; }

int main(){
  // ====== 可配置参数 ======
  const int n = 3;                 // 参与方数量 P1,P2,P3
  const int k = 2;                 // 阈值（示例设置 k=2）
  const double eps_okvs = 0.07;    // OKVS 负载
  const size_t w = 64;             // RB-OKVS 带宽（占位）
  const double eps_hash = 1.3;     // 哈希表负载
  const u64 salt_tag = 0x13572468abcdef99ULL;  // 标签盐
  const u64 seed_pos = 0x24681357abcdef01ULL;  // 位置哈希种子

  // ====== 构造每方集合 X_i（示例：部分重叠） ======
  std::vector<std::vector<Block128>> Xs(n+1);
  // P1
  Xs[1] = { X(1,10), X(2,20), X(3,30), X(4,40) };
  // P2
  Xs[2] = { X(2,20), X(3,30), X(5,50) };
  // P3
  Xs[3] = { X(3,30), X(6,60) };

  // ====== S12: 每方构造 (x, f_x(i)) 与 tag ======
  std::vector<std::vector<KV>> kv_all(n+1);
  std::vector<std::vector<Tag128>> tag_all(n+1);
  for(int i=1;i<=n;++i){
    auto alpha_i = hash_to_field_u64(i);
    for(auto& x : Xs[i]){
      PRG prg(mix_seed(x, /*salt=*/0xC0FFEEULL ^ i)); // 为每个 x 派生子种子
      auto fxi = poly_eval_fx_at_i(x, alpha_i, k, prg);
      kv_all[i].push_back( KV{x, fxi} );
      tag_all[i].push_back( tag_of(x, salt_tag) );
    }
  }

  // ====== S13: OKVS 编码 & “交换”（本示例直接内存中共享） ======
  std::vector<RBOKVS> okvs(n+1);
  std::vector<size_t> ni(n+1);
  for(int i=1;i<=n;++i){
    ni[i] = Xs[i].size();
    OKVSParams p;
    p.m = size_t((1.0+eps_okvs)*ni[i]) + 1;
    p.w = w; p.seed_r1 = 1234+i; p.seed_r2 = 5678+i;
    okvs[i] = RBOKVS::Encode(kv_all[i], p);
  }

  // ====== S14: 解码 σ_{i,j}^γ ======
  // 构造 T_i（哈希表），B = eps_hash * M
  size_t M = 0; for(int i=1;i<=n;++i) M = std::max(M, ni[i]);
  size_t B = size_t(eps_hash * M) + 1;
  std::vector<HashTableTi> Ts(n+1, HashTableTi{std::vector<Bucket>(B)});

  for(int i=1;i<=n;++i){
    for(size_t j=0;j<ni[i];++j){
      const auto& x = kv_all[i][j].key;
      const auto& fxi = kv_all[i][j].val;
      std::vector<std::pair<int, Block128>> sigmas;
      for(int g=1; g<=n; ++g) if(g!=i){
        auto sig = okvs[g].Decode(x);
        sigmas.emplace_back(g, sig);
      }
      insert_element_Ti(Ts[i], B, n, i, x, tag_all[i][j], fxi, sigmas, seed_pos);
    }
  }

  // ====== S31–S33: 设 P_n 为终端聚合方，聚合同桶份额，判断交集 ======
  // 简化的“验证”：若同一 tag 在 n 个 T_i 中累计命中份额数 >= k（且包含至少一个“本方真值”），认为命中交集
  // 注意：这里为了演示，直接以 okvs[?] 的 Decode 来模拟“真假份额”。真实实现需插值验证。
  std::set<u64> result_hash;             // 用 lo 作为可视化键（演示）
  std::vector<Block128> result_blocks;   // 收集命中元素（演示按 P_n 自己的集合对齐）

  int agg = n; // 终端聚合方 = P_n
  for(size_t eta=0; eta<B; ++eta){
    // 收集所有方该桶的条目
    std::vector<Share> pool;
    for(int i=1;i<=n;++i){
      for(auto &sh : Ts[i].table[eta].items) pool.push_back(sh);
    }
    // 按 tag 聚类
    std::unordered_map<u64, std::vector<Share>> by_tag;
    auto tag_key = [](const Tag128& t)->u64 { return t.hi ^ (t.lo<<1); };
    for(auto &sh: pool) by_tag[tag_key(sh.tag)].push_back(sh);

    // 计数并判定
    for(auto& kv : by_tag){
      auto &vec = kv.second;
      // 去重（同一 party 取一份）
      std::set<int> parties;
      for(auto &sh: vec) parties.insert(sh.party_id);
      if((int)parties.size() >= k){
        result_hash.insert(kv.first);
      }
    }
  }

  // 为了可视化：从 P_n 的集合中找出这些 tag 对应的 x（演示）
  std::unordered_map<u64, Block128> tag2x;
  for(size_t j=0;j<ni[agg];++j){
    auto t = tag_all[agg][j];
    u64 key = (t.hi ^ (t.lo<<1));
    tag2x[key] = kv_all[agg][j].key;
  }
  for(auto h : result_hash){
    auto it = tag2x.find(h);
    if(it != tag2x.end()) result_blocks.push_back(it->second);
  }

  // ====== 输出到终端 ======
  std::cout << "=== Estimated intersection (by tag threshold >= " << k << ") ===\n";
  for(size_t i=0;i<result_blocks.size();++i){
    std::cout << i << ": (" << result_blocks[i].hi << "," << result_blocks[i].lo << ")\n";
  }
  std::cout << "Total: " << result_blocks.size() << "\n";

  // ====== 写 CSV 以便画图 ======
  write_csv("result.csv", result_blocks);
  std::cout << "CSV written to result.csv\n";
  return 0;
}
