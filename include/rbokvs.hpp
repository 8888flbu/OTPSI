#pragma once
#include "types.hpp"
#include <unordered_map>

// 占位 RB-OKVS：用 unordered_map 代替（可运行、便于先打通流程）
// 之后你把 Encode/Decode 改成真正的 RB-OKVS 即可。
struct RBOKVS {
  OKVSParams p;
  // 演示数据：直接存 k->v 的表（真实实现会是编码向量 s）
  std::unordered_map<u64, Block128> kv;

  static RBOKVS Encode(const std::vector<KV>& kvs, const OKVSParams& p);
  Block128 Decode(const Block128& key) const;
};

// 简易 key 哈希（仅演示；真实请用 BLAKE3）
inline u64 key_hash64(const Block128& k){
  return k.hi ^ (k.lo*0x9e3779b97f4a7c15ULL);
}
