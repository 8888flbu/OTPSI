#pragma once
#include "types.hpp"   // 已经有 OKVSParams、Block128、KV
#include <vector>

// RB-OKVS 存储结构
struct RBOKVS {
    OKVSParams p;
    std::vector<Block128> S;  // 长度 m 的存储向量

    // 编码：把若干 (key, value) 映射到 S
    static RBOKVS Encode(const std::vector<KV>& kvs, const OKVSParams& p);

    // 解码：从 key 恢复 value
    Block128 Decode(const Block128& key) const;
};

// ==== 哈希辅助函数 ====

// 将 key 哈希到 [0, m-w]
inline size_t H1(const Block128& k, const OKVSParams& p) {
    uint64_t h = k.hi ^ (k.lo * 0x9e3779b97f4a7c15ULL ^ p.seed_r1);
    return static_cast<size_t>(h % (p.m - p.w + 1));
}

// 将 key 映射为 w 比特的布尔向量 u
inline std::vector<uint8_t> H2(const Block128& k, const OKVSParams& p) {
    uint64_t h = k.lo ^ (k.hi * 0x9e3779b97f4a7c15ULL ^ p.seed_r2);
    std::vector<uint8_t> u(p.w);
    for (uint32_t j = 0; j < p.w; ++j) {
        u[j] = (h >> (j % 64)) & 1; // 简化：从 64bit 循环取位
    }
    return u;
}

