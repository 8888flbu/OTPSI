#pragma once
#include "types.hpp"
#include <fstream>

// 非网络场景：简单把结果写 CSV 便于可视化
inline void write_csv(const std::string& path, const std::vector<Block128>& xs){
  std::ofstream ofs(path);
  ofs << "idx,hi,lo\n";
  for(size_t i=0;i<xs.size();++i) ofs << i << "," << xs[i].hi << "," << xs[i].lo << "\n";
}
