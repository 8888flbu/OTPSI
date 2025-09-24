#pragma once
#include "types.hpp"

// 占位实现：add = XOR；mul = 仅演示用途（非真 GF 乘法！）
inline Block128 gf_add(const Block128& a, const Block128& b){
  return Block128{a.hi ^ b.hi, a.lo ^ b.lo};
}

// TODO: 替换为真正的 GF(2^128) 乘法（CLMUL/PMULL 最佳）
inline Block128 gf_mul(const Block128& a, const Block128& b){
  // 演示：用 64-bit 乘法简化（不是正式安全实现）
  u64 lo = (a.lo + 0x9e3779b97f4a7c15ULL) * (b.lo | 1ULL);
  u64 hi = (a.hi + 0x517cc1b727220a95ULL) * (b.hi | 1ULL);
  return Block128{hi, lo};
}

// 把整数 i 映射到“域元素”
inline Block128 hash_to_field_u64(u64 i){
  // 演示：SplitMix64 双次（建议你改为 BLAKE3 派生 128-bit）
  auto mix = [](u64 x){
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
  };
  return Block128{mix(i), mix(i ^ 0xD1B54A32D192ED03ULL)};
}
