#include "rbokvs.hpp"
#include <algorithm>
#include <unordered_map>
#include <cstdint>
#include <chrono>
#include <fstream>

void benchmark_okvs(const std::vector<KV>& kvs, const OKVSParams& p) {
    using namespace std::chrono;

    auto start = high_resolution_clock::now();
    auto okvs = RBOKVS::Encode(kvs, p);
    auto end = high_resolution_clock::now();
    double ms = duration<double, std::milli>(end - start).count();

    std::ofstream fout("okvs_bench.csv", std::ios::app);
    fout << "encode," << kvs.size() << "," << p.m << "," 
         << (double)kvs.size() / p.m << "," << ms << "\n";
}


// ============ 小工具：GF(2^128) 加法 = 128bit XOR ============
static inline void xor_inplace(Block128& a, const Block128& b) {
    a.hi ^= b.hi; a.lo ^= b.lo;
}
static inline bool is_zero128(const Block128& x) {
    return (x.hi | x.lo) == 0ull;
}

// ============ 极小概率“全 0 带”兜底 ============
static inline void ensure_nonzero(std::vector<uint8_t>& u) {
    bool all0 = true;
    for (auto b : u) { if (b) { all0 = false; break; } }
    if (all0 && !u.empty()) u[0] = 1;
}

// ============ 简单可复现 PRG：splitmix64 ============
static inline uint64_t splitmix64(uint64_t& x) {
    uint64_t z = (x += 0x9e3779b97f4a7c15ull);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
}
static inline Block128 prg_block(uint64_t k1, uint64_t k2) {
    uint64_t s = k1 ^ (k2 * 0x9e3779b97f4a7c15ull);
    Block128 b{};
    b.hi = splitmix64(s);
    b.lo = splitmix64(s);
    return b;
}

// ============ 行结构 ============
struct RowBand {
    size_t a;
    std::vector<uint8_t> u;
    Block128 v;

    size_t first_one() const {
        for (size_t j = 0; j < u.size(); ++j) if (u[j]) return j;
        return u.size();
    }
    void xor_with(const RowBand& other) {
        const size_t W = u.size();
        for (size_t j = 0; j < W; ++j) u[j] ^= other.u[j];
        xor_inplace(v, other.v);
    }
};

// ============ 编码 ============
RBOKVS RBOKVS::Encode(const std::vector<KV>& kvs, const OKVSParams& p) {
    using clk = std::chrono::high_resolution_clock;
    auto t0 = clk::now();

    RBOKVS out; out.p = p;
    const size_t m = p.m;
    const uint32_t w = p.w;

    if (m == 0 || w == 0 || m <= w) {
        out.S.assign(std::max<size_t>(1, m ? m : 1), Block128{0,0});
        for (size_t i = 0; i < out.S.size(); ++i)
            out.S[i] = prg_block(p.seed_r1 ^ (uint64_t)i, p.seed_r2 + (uint64_t)m);
        return out;
    }

    // 1) 构造行
    std::vector<RowBand> rows; rows.reserve(kvs.size());
    for (const auto& e : kvs) {
        RowBand r;
        r.a = H1(e.key, p);
        r.u = H2(e.key, p);
        ensure_nonzero(r.u);
        r.v = e.val;
        rows.emplace_back(std::move(r));
    }

    // 2) 排序
    auto lead_col = [&](const RowBand& r)->size_t {
        size_t j = r.first_one();
        return (j < r.u.size()) ? (r.a + j) : (size_t)-1;
    };
    std::sort(rows.begin(), rows.end(), [&](const RowBand& x, const RowBand& y){
        return lead_col(x) < lead_col(y);
    });

    // 3) 消元
    std::vector<RowBand> basis; basis.reserve(rows.size());
    std::unordered_map<size_t, size_t> pivot_at;
    pivot_at.reserve(rows.size() * 2);

    for (auto& r : rows) {
        for (size_t off = 0; off < w; ++off) {
            if (!r.u[off]) continue;
            const size_t col = r.a + off;
            auto it = pivot_at.find(col);
            if (it != pivot_at.end()) {
                r.xor_with(basis[it->second]);
            }
        }
        size_t j = r.first_one();
        if (j == w) {
            if (!is_zero128(r.v)) {
                out.S.assign(m, Block128{0,0});
                for (size_t i = 0; i < m; ++i)
                    out.S[i] = prg_block(p.seed_r1 + (uint64_t)i, p.seed_r2 ^ 0xA5A5A5A5A5A5A5A5ull);
                return out;
            }
            continue;
        }
        const size_t col = r.a + j;
        pivot_at[col] = basis.size();
        basis.emplace_back(std::move(r));
    }

    // 4) 自由列随机化
    out.S.assign(m, Block128{0,0});
    std::vector<uint8_t> is_free(m, 1);
    for (auto& pr : pivot_at) is_free[pr.first] = 0;
    for (size_t col = 0; col < m; ++col) {
        if (is_free[col]) {
            out.S[col] = prg_block(p.seed_r1 ^ (uint64_t)(0x1111111111111111ull + col),
                                   p.seed_r2 ^ (uint64_t)(0x2222222222222222ull + col));
        }
    }

    // 5) 回代
    std::vector<size_t> piv_cols; piv_cols.reserve(pivot_at.size());
    for (auto& pr : pivot_at) piv_cols.push_back(pr.first);
    std::sort(piv_cols.begin(), piv_cols.end(), std::greater<size_t>());

    for (size_t pc : piv_cols) {
        const RowBand& r = basis[pivot_at[pc]];
        size_t jstar = r.first_one();
        Block128 acc = r.v;
        for (size_t j = 0; j < w; ++j) {
            if (!r.u[j] || j == jstar) continue;
            const size_t col = r.a + j;
            xor_inplace(acc, out.S[col]);
        }
        out.S[pc] = acc;
    }

    auto t1 = clk::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // 写入 CSV（追加模式）
    std::ofstream ofs("okvs_bench.csv", std::ios::app);
    ofs << "encode," << kvs.size() << "," << m << "," << (double)kvs.size()/m
        << "," << ms << "\n";
    ofs.close();

    return out;
}

// ============ 解码 ============
Block128 RBOKVS::Decode(const Block128& key) const {
    using clk = std::chrono::high_resolution_clock;
    auto t0 = clk::now();

    const size_t a = H1(key, p);
    auto u = H2(key, p);
    ensure_nonzero(u);
    Block128 acc{0,0};
    for (uint32_t j = 0; j < p.w; ++j) {
        if (u[j]) xor_inplace(acc, S[a + j]);
    }

    auto t1 = clk::now();
    double us = std::chrono::duration<double, std::micro>(t1 - t0).count();

    // 追加写 CSV
    std::ofstream ofs("okvs_bench.csv", std::ios::app);
    ofs << "decode,1," << p.m << "," << (double)1/p.m << "," << us/1000.0 << "\n";
    ofs.close();

    return acc;
}

