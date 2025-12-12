#pragma once
#include <cstdint>
#include "types.hpp"

// 假定 types.hpp 定义：using u64 = uint64_t; struct Block128{ u64 hi, lo; };

namespace gf128 {

// ====== 基本工具 ======
static inline Block128 zero() { return Block128{0,0}; }
static inline Block128 one()  { return Block128{0,1}; }

static inline bool is_zero(const Block128& a){
    return (a.hi | a.lo) == 0ULL;
}

// 加法 = XOR
static inline Block128 add(const Block128& a, const Block128& b){
    return Block128{ a.hi ^ b.hi, a.lo ^ b.lo };
}
static inline Block128 sub(const Block128& a, const Block128& b){
    return add(a,b); // 同加法
}

// ====== 进位无关 64x64 乘法：可移植慢路径 ======
static inline void clmul64(u64 a, u64 b, u64& hi, u64& lo){
    unsigned __int128 acc = 0;
    for (int i=0;i<64;i++){
        if ((b >> i) & 1ULL) acc ^= ( (unsigned __int128)a << i );
    }
    lo = (u64)acc;
    hi = (u64)(acc >> 64);
}

// 128x128 的 carry-less 乘法（Karatsuba 合成 256bit）
static inline void clmul128(const Block128& A, const Block128& B,
                            Block128& lo128, Block128& hi128){
    u64 a0 = A.lo, a1 = A.hi;
    u64 b0 = B.lo, b1 = B.hi;

    u64 p00_lo, p00_hi, p11_lo, p11_hi, pm_lo, pm_hi;
    clmul64(a0, b0, p00_hi, p00_lo);             // p00 = a0*b0
    clmul64(a1, b1, p11_hi, p11_lo);             // p11 = a1*b1
    u64 t0 = a0 ^ a1, t1 = b0 ^ b1;
    u64 pm_lo_lo, pm_lo_hi;
    clmul64(t0, t1, pm_lo_hi, pm_lo_lo);         // pm = (a0^a1)*(b0^b1)

    // mid = pm ^ p00 ^ p11
    u64 mid_lo = pm_lo_lo ^ p00_lo ^ p11_lo;
    u64 mid_hi = pm_lo_hi ^ p00_hi ^ p11_hi;

    // 组装 256bit: r3:r2:r1:r0
    u64 r0 = p00_lo;
    u64 r1 = p00_hi ^ mid_lo;
    u64 r2 = p11_lo ^ mid_hi;
    u64 r3 = p11_hi;

    lo128 = Block128{ r1, r0 };   // 低 128 位
    hi128 = Block128{ r3, r2 };   // 高 128 位
}

// ====== 模多项式约简：x^128 + x^7 + x^2 + x + 1（GHASH） ======
static inline Block128 reduce_256(const Block128& hi, const Block128& lo){
    Block128 res = lo;

    auto xor_bit = [](Block128& v, int pos){
        if (pos < 0 || pos >= 128) return;
        if (pos < 64) v.lo ^= (1ULL << pos);
        else          v.hi ^= (1ULL << (pos - 64));
    };

    auto fold_word = [&](u64 word, int base){
        while (word){
            int t = __builtin_ctzll(word); // 取最低 1 的位置
            int k = base + t;              // hi 的第 k 位为 1
            xor_bit(res, k);
            xor_bit(res, k+1);
            xor_bit(res, k+2);
            xor_bit(res, k+7);
            word &= (word - 1);
        }
    };

    fold_word(hi.lo, 0);
    fold_word(hi.hi, 64);
    return res;
}

// 乘法：clmul128 → 约简
static inline Block128 mul(const Block128& a, const Block128& b){
    if (is_zero(a) || is_zero(b)) return zero();
    Block128 lo128{}, hi128{};
    clmul128(a,b,lo128,hi128);
    return reduce_256(hi128, lo128);
}

// 平方（可直接 mul(a,a)）
static inline Block128 square(const Block128& a){
    return mul(a,a);
}

// 求逆：a^(2^128 - 2)
static inline Block128 inv(const Block128& a){
    if (is_zero(a)) return zero(); // 0 无逆
    Block128 result = one();
    Block128 base   = a;
    // 处理位 127..1 为 1，位 0 为 0
    for (int i = 127; i >= 1; --i){
        result = mul(result, base);
        base   = square(base);
    }
    base = square(base); // 推进最低位（0）的一次平方
    return result;
}

// 将 u64 编码到 GF(2^128)（低 64 位放值，高位 0）
static inline Block128 from_u64(u64 x){ return Block128{0, x}; }

// 把整数 i 映射到 GF(2^128) 元素（稳定可复现哈希）
inline Block128 hash_to_field_u64(u64 i){
  auto mix = [](u64 x){
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
  };
  return Block128{ mix(i), mix(i ^ 0xD1B54A32D192ED03ULL) };
}

} // namespace gf128

// ====== 兼容层（为旧代码保留全局函数名） ======
static inline Block128 gf_add(const Block128& a, const Block128& b) { return gf128::add(a,b); }
static inline Block128 gf_mul(const Block128& a, const Block128& b) { return gf128::mul(a,b); }
static inline Block128 gf_square(const Block128& a)                 { return gf128::square(a); }
static inline Block128 gf_inv(const Block128& a)                    { return gf128::inv(a); }
static inline Block128 from_u64(u64 x)                              { return gf128::from_u64(x); }
static inline Block128 hash_to_field_u64(u64 i)                     { return gf128::hash_to_field_u64(i); }

