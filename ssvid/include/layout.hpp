#pragma once
#include "types.hpp"
#include <algorithm>

inline u64 h_pos(const Block128& x, u64 dom, u64 seed, u64 B){
  // 简化位置哈希（真实请用 BLAKE3 域分离）
  u64 h = (x.hi ^ (x.lo<<1) ^ dom ^ seed) * 0x9e3779b97f4a7c15ULL;
  h ^= (h>>33); h *= 0xff51afd7ed558ccdULL; h ^= (h>>33);
  return h % B;
}

inline std::vector<size_t> positions_sorted_n(const Block128& x, int n, u64 seed_pos, size_t B){
  std::vector<size_t> I; I.reserve(n);
  for(int ell=1; ell<=n; ++ell) I.push_back(h_pos(x, ell, seed_pos, B));
  std::sort(I.begin(), I.end());
  return I;
}

// 插入规则（与你的 S23 对齐）
inline void insert_element_Ti(
  HashTableTi& Ti, size_t B, int n, int i,
  const Block128& x, const Tag128& tag_x,
  const Block128& fxi,
  const std::vector<std::pair<int, Block128>>& sigmas_gamma,
  u64 seed_pos
){
  auto I = positions_sorted_n(x, n, seed_pos, B);
  size_t idx_self = I[( (i + i - 1) % n )];
  Ti.table[idx_self].items.push_back( Share{i, fxi, tag_x} );
  for(auto& [g, sig] : sigmas_gamma){
    size_t idx_g = I[( (g + i - 1) % n )];
    Ti.table[idx_g].items.push_back( Share{g, sig, tag_x} );
  }
}
