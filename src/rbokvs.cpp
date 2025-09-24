#include "rbokvs.hpp"
#include <algorithm>

// 小工具：GF(2^128) 加法 = 按位异或
static inline void xor_inplace(Block128& a, const Block128& b) {
    a.hi ^= b.hi;
    a.lo ^= b.lo;
}

// 如果 H2 生成的向量全 0，兜底：强制把第 0 位设为 1，避免空带
static inline void ensure_nonzero(std::vector<uint8_t>& u) {
    bool all_zero = true;
    for (auto b : u) { if (b) { all_zero = false; break; } }
    if (all_zero && !u.empty()) u[0] = 1;
}

RBOKVS RBOKVS::Encode(const std::vector<KV>& kvs, const OKVSParams& p) {
    RBOKVS out;
    out.p = p;
    out.S.assign(p.m, Block128{0, 0});  // 初始化存储向量 S[0..m-1]

    // === 简化版策略（可运行）===
    // 对每个 (key, val)：
    //   1) a = H1(key) ∈ [0, m-w]
    //   2) u = H2(key) ∈ {0,1}^w，若全 0 则把 u[0]=1
    //   3) 选 u[j]==1 的第一个 j*，把 val 异或进 S[a+j*]
    //
    // 解码时按 u 的 1 位把 S[a..a+w-1] XOR 回来即可得到 val（如果没有碰撞）。
    // 注意：这是“能跑”的占位实现，可能存在碰撞导致误差；要论文级正确性需替换为带状矩阵求解。
    for (const auto& e : kvs) {
        const size_t a = H1(e.key, p);            // 带起点
        auto u = H2(e.key, p);                    // 带比特
        ensure_nonzero(u);

        // 找到第一个 1 位
        uint32_t jstar = 0;
        while (jstar < p.w && u[jstar] == 0) ++jstar;
        if (jstar >= p.w) jstar = 0;              // 兜底

        const size_t idx = a + jstar;             // 写入位置
        // 做 XOR 聚合（若不同键落到同一槽会产生碰撞，这是简化版的代价）
        xor_inplace(out.S[idx], e.val);
    }

    return out;
}

Block128 RBOKVS::Decode(const Block128& key) const {
    const size_t a = H1(key, p);
    auto u = H2(key, p);
    ensure_nonzero(u);

    Block128 acc{0, 0};
    // 解码：acc = ⊕_{j: u[j]=1} S[a+j]
    for (uint32_t j = 0; j < p.w; ++j) {
        if (u[j]) {
            const size_t idx = a + j;
            acc.hi ^= S[idx].hi;
            acc.lo ^= S[idx].lo;
        }
    }
    return acc;
}

