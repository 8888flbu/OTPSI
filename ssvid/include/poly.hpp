#pragma once
#include "types.hpp"
#include "gf128.hpp"
#include "hash_prg.hpp"

// f_x(t) = x ⊕ r1 t ⊕ ... ⊕ r_{k-1} t^{k-1}，Horner 计算 f_x(i)
inline Block128 poly_eval_fx_at_i(const Block128& x, const Block128& alpha_i, int k, PRG& prg_for_x){
  if(k<=1) return x;
  std::vector<Block128> r(k-1);
  for(int d=0; d<k-1; ++d) r[d] = prg_for_x.next_block128();
  Block128 acc = r.back();
  for(int d=k-3; d>=0; --d) acc = gf_add(gf_mul(acc, alpha_i), r[d]);
  return gf_add(gf_mul(acc, alpha_i), x);
}
