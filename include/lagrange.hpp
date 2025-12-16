#pragma once
#include <vector>
#include "types.hpp"
#include "gf128.hpp"

namespace gf128 {

// 拉格朗日插值：在 GF(2^128) 上，用 k 个点 (x_i, y_i) 计算 f(x0)
static inline Block128 lagrange_at_gf128(
    const std::vector<std::pair<Block128, Block128>>& pts,
    const Block128& x0)
{
    using namespace gf128;
    Block128 acc = zero();

    for (size_t i = 0; i < pts.size(); ++i){
        const Block128 xi = pts[i].first;
        const Block128 yi = pts[i].second;

        // 计算 L_i(x0) = Π_{j≠i} (x0 - x_j)/(x_i - x_j)
        // 在 GF(2^m) 中 - == +，因此 (x0 - xj) == (x0 ^ xj)
        Block128 num = one();
        Block128 den = one();
        for (size_t j = 0; j < pts.size(); ++j){
            if (j == i) continue;
            const Block128 xj = pts[j].first;
            Block128 t_num = add(x0, xj);
            Block128 t_den = add(xi, xj);
            num = mul(num, t_num);
            den = mul(den, t_den);
        }
        Block128 Li = mul(num, inv(den));
        acc = add(acc, mul(yi, Li));
    }
    return acc; // f(x0)
}

// 便捷：把 u64 作为横坐标（参与方编号）→ GF(2^128)
static inline Block128 x_from_party_id_u64(u64 party_id){
    return from_u64(party_id);
}

} // namespace gf128
