#include "rbokvs.hpp"

RBOKVS RBOKVS::Encode(const std::vector<KV>& kvs, const OKVSParams& p){
  RBOKVS out; out.p = p;
  for(auto& e: kvs) out.kv[key_hash64(e.key)] = e.val;
  return out;
}
Block128 RBOKVS::Decode(const Block128& key) const{
  auto it = kv.find(key_hash64(key));
  if(it!=kv.end()) return it->second;
  // 未命中返回“伪随机”块（演示：用 hash 返回）
  u64 h = key_hash64(key);
  return Block128{ h ^ 0xA5A5A5A5A5A5A5A5ULL, h ^ 0x5A5A5A5A5A5A5A5AULL };
}
