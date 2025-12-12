cmake_minimum_required(VERSION 3.20)

# 同时启用 C / C++ / ASM（BLAKE3 是 C 源码）
project(ssvid C CXX ASM)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ===================== BLAKE3 路径探测 =====================
set(BLAKE3_DIR ${CMAKE_SOURCE_DIR}/third_party/blake3)
if(EXISTS ${BLAKE3_DIR}/c/blake3.h)
  set(BLAKE3_SRC_DIR ${BLAKE3_DIR}/c)   # 官方 C 参考实现目录
else()
  set(BLAKE3_SRC_DIR ${BLAKE3_DIR})     # 兼容旧布局
endif()
# ...前略...

# （此前已设置 BLAKE3_SRC_DIR 指到 third_party/blake3/c）
add_library(core
  src/gf128.cpp
  src/hash_prg.cpp
  src/poly.cpp
  src/rbokvs.cpp
  src/layout.cpp
  src/wire.cpp
  ${BLAKE3_SRC_DIR}/blake3.c
  ${BLAKE3_SRC_DIR}/blake3_dispatch.c
  ${BLAKE3_SRC_DIR}/blake3_portable.c
)

target_include_directories(core
  PUBLIC  ${CMAKE_SOURCE_DIR}/include
          ${BLAKE3_SRC_DIR}
)

if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
  enable_language(ASM)
  foreach(asmfile
      blake3_sse2_x86-64_unix.S
      blake3_sse41_x86-64_unix.S   # ← 新增这一行
      blake3_avx2_x86-64_unix.S
      blake3_avx512_x86-64_unix.S)
    if(EXISTS ${BLAKE3_SRC_DIR}/${asmfile})
      target_sources(core PRIVATE ${BLAKE3_SRC_DIR}/${asmfile})
    endif()
  endforeach()
endif()

add_executable(party src/main_party.cpp)
target_link_libraries(party PRIVATE core)



