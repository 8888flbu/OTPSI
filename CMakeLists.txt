cmake_minimum_required(VERSION 3.20)
project(ssvid CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_library(core
  src/gf128.cpp
  src/hash_prg.cpp
  src/poly.cpp
  src/rbokvs.cpp
  src/layout.cpp
  src/wire.cpp
)
target_include_directories(core PUBLIC include)
add_executable(party src/main_party.cpp)
target_link_libraries(party PRIVATE core)
