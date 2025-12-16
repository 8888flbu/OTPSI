#pragma once
#include "types.hpp"
#include <random>
#include <cstring>

// 占位 PRG：std::mt19937_64（可重复）；建议改为 BLAKE3-XOF
struct PRG {
  std::mt19937_64 eng;
  explicit PRG(u64 seed){ eng.seed(seed); }
  Block128 next_block128(){
    u64 a = eng(); u64 b = eng();
    return Block128{a, b};
  }
};

inline u64 mix_seed(const Block128& blk, u64 salt){
  return blk.hi ^ (blk.lo<<1) ^ (salt * 0x9e3779b97f4a7c15ULL);
}

// 演示标签：把 x 与盐混合后映射到 Tag128（建议改 BLAKE3）
inline Tag128 tag_of(const Block128& x, u64 salt){
  PRG prg(mix_seed(x, salt));
  auto b = prg.next_block128();
  return Tag128{b.hi, b.lo};
}
