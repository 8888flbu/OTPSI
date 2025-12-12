#pragma once
#include <array>
#include <cstdint>
#include <vector>
#include <string>
#include <utility>

using u8 = uint8_t; using u64 = uint64_t;

struct Block128 { u64 hi{0}, lo{0}; };         // 128-bit 块（GF(2^128) 占位）
struct Tag128   { u64 hi{0}, lo{0}; };         // 128-bit 标签（哈希到此）

inline bool operator==(const Tag128& a, const Tag128& b){ return a.hi==b.hi && a.lo==b.lo; }

struct KV { Block128 key; Block128 val; };     // (x, f_x(i))
struct Share { int party_id; Block128 fx_i; Tag128 tag; };

struct OKVSParams { size_t m{0}; size_t w{64}; u64 seed_r1{0}; u64 seed_r2{0}; };

struct Bucket { std::vector<Share> items; };
struct HashTableTi { std::vector<Bucket> table; /* size B */ };

// 便捷工具
inline Block128 make_block(u64 hi, u64 lo){ return Block128{hi, lo}; }
inline Tag128   make_tag  (u64 hi, u64 lo){ return Tag128  {hi, lo}; }
