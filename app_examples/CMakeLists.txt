cmake_minimum_required(VERSION 3.15)

# The SDK (setup_sdk.cmake) configures effcc as the compiler. Declare no languages here
# so CMake doesn't probe for a host compiler first: on Windows that probe picks up MSVC or
# MinGW and injects host link flags (-luuid, --out-implib, ...) that the effcc linker
# rejects. CMAKE_SYSTEM_NAME must be set before project() for the same reason.
set(CMAKE_SYSTEM_NAME Generic)
project(app_examples LANGUAGES NONE)

# Uncomment to enable hand-optimized app versions:
#
# set(EFF_BLD_HAND_OPTIMIZED ON)

# EVK Requires UART peripheral 3
set(EFF_STDIO_PORT 3)

# Pass the EFF SDK location with -DEFF_SDK_ROOT_DIR=<path>
if(NOT DEFINED EFF_SDK_ROOT_DIR)
  message(
    FATAL_ERROR
      "Set EFF_SDK_ROOT_DIR to the EFF SDK, for example -DEFF_SDK_ROOT_DIR=<effcc package>/sdk"
  )
endif()
file(REAL_PATH "${EFF_SDK_ROOT_DIR}" EFF_SDK_ROOT_DIR)
include(${EFF_SDK_ROOT_DIR}/setup_sdk.cmake)
add_subdirectory(${EFF_SDK_ROOT_DIR} sdk)

if(EFF_BLD_HAND_OPTIMIZED)
  add_compile_definitions(EFF_BLD_HAND_OPTIMIZED)
endif()

add_subdirectory(bfs)
add_subdirectory(biquad_filter)
add_subdirectory(cholesky_decomp)
add_subdirectory(conv3x3)
add_subdirectory(conv3x3_dw)
add_subdirectory(conv3x3xn)
add_subdirectory(conv5x5)
add_subdirectory(crc32)
add_subdirectory(dijkstras)
add_subdirectory(dmm)
add_subdirectory(dmm_i8)
add_subdirectory(dmv)
add_subdirectory(dmv_fp)
add_subdirectory(fft)
add_subdirectory(fir)
add_subdirectory(hello)
add_subdirectory(isqrt)
add_subdirectory(jpeg)
add_subdirectory(knapsack)
add_subdirectory(ldpc)
add_subdirectory(quickstart)
add_subdirectory(levenshtein)
add_subdirectory(qr_decomp)
add_subdirectory(recursive_factorial)
add_subdirectory(smv)
add_subdirectory(sort)
add_subdirectory(spadd)
add_subdirectory(sparse)
add_subdirectory(sqrt_newton)
