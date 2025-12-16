#pragma once
#include "types.hpp"   // 已经有 OKVSParams、Block128、KV
#include <vector>
#include <cstdint>
#include <cstring>

// 引入 BLAKE3 头（C 接口）
extern "C" {
#include "blake3.h"
}

// ================= RB-OKVS 存储结构 =================
struct RBOKVS {
    OKVSParams p;
    std::vector<Block128> S;  // 长度 m 的存储向量

// *** COMM *** 新增：计算 OKVS 的字节大小（用于通信量统计）
    size_t byte_size() const {
        return S.size() * sizeof(Block128);
    }

    // 编码：把若干 (key, value) 映射到 S
    static RBOKVS Encode(const std::vector<KV>& kvs, const OKVSParams& p);
    // 解码：从 key 恢复 value
    Block128 Decode(const Block128& key) const;
};

// ------------- 工具：Block128 序列化（大端） -------------
inline void ser_block128_be(const Block128& k, uint8_t out[16]) {
    auto wr64be = [](uint64_t x, uint8_t* b){
        b[0]=(uint8_t)(x>>56); b[1]=(uint8_t)(x>>48); b[2]=(uint8_t)(x>>40); b[3]=(uint8_t)(x>>32);
        b[4]=(uint8_t)(x>>24); b[5]=(uint8_t)(x>>16); b[6]=(uint8_t)(x>>8);  b[7]=(uint8_t)(x);
    };
    wr64be(k.hi, out);
    wr64be(k.lo, out+8);
}

// ------------- 用 64-bit 种子构造 32B BLAKE3 密钥 -------------
inline void fill_key_from_seed(uint64_t seed, uint8_t key[32]) {
    uint64_t s0 = seed ^ 0x9e3779b97f4a7c15ULL;
    uint64_t s1 = seed ^ 0xbf58476d1ce4e5b9ULL;
    uint64_t s2 = seed ^ 0x94d049bb133111ebULL;
    uint64_t s3 = seed ^ 0x2545F4914F6CDD1DULL;
    std::memcpy(key +  0, &s0, 8);
    std::memcpy(key +  8, &s1, 8);
    std::memcpy(key + 16, &s2, 8);
    std::memcpy(key + 24, &s3, 8);
}

// ================= 安全版 H1/H2（BLAKE3 keyed） =================

// H1：把 key 均匀映射到 [0, m-w] 区间
inline size_t H1(const Block128& k, const OKVSParams& p) {
    const size_t range = (p.m > p.w) ? (p.m - p.w + 1) : 1;

    uint8_t key32[32]; fill_key_from_seed(p.seed_r1, key32);
    uint8_t in[16];    ser_block128_be(k, in);
    uint8_t out8[8];

    blake3_hasher hasher;
    blake3_hasher_init_keyed(&hasher, key32);
    blake3_hasher_update(&hasher, in, 16);
    blake3_hasher_finalize_seek(&hasher, 0, out8, sizeof(out8));

    // 按 LE 读 64bit（只要 Encode/Decode 一致即可）
    uint64_t x = ((uint64_t)out8[0])       |
                 ((uint64_t)out8[1] << 8)  |
                 ((uint64_t)out8[2] << 16) |
                 ((uint64_t)out8[3] << 24) |
                 ((uint64_t)out8[4] << 32) |
                 ((uint64_t)out8[5] << 40) |
                 ((uint64_t)out8[6] << 48) |
                 ((uint64_t)out8[7] << 56);

    return (size_t)(x % range);
}

// H2：把 key 映射为长度为 w 的 {0,1} 向量（任意 w，XOF 流按位取）
inline std::vector<uint8_t> H2(const Block128& k, const OKVSParams& p) {
    const uint32_t w = p.w;
    std::vector<uint8_t> bits(w);

    uint8_t key32[32]; fill_key_from_seed(p.seed_r2, key32);
    uint8_t in[16];    ser_block128_be(k, in);

    blake3_hasher hasher;
    blake3_hasher_init_keyed(&hasher, key32);
    blake3_hasher_update(&hasher, in, 16);

    const size_t nbytes = (w + 7) / 8;
    std::vector<uint8_t> stream(nbytes, 0);
    blake3_hasher_finalize_seek(&hasher, 0, stream.data(), stream.size());

    // 展开到位（第 j 位取 stream[j>>3] 的 (j&7) 位）
    for (uint32_t j = 0; j < w; ++j) {
        bits[j] = (stream[j >> 3] >> (j & 7)) & 1;
    }

    // 极小概率全 0（≈2^-w），做兜底以避免退化行
    if (w > 0) {
        bool all_zero = true;
        for (auto b : bits) { if (b) { all_zero = false; break; } }
        if (all_zero) bits[0] = 1;
    }
    return bits;
}

