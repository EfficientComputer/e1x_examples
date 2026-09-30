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

# The EFF SDK is the sdk folder of the installed effcc package. Pass it with
# -DEFF_SDK_ROOT_DIR=<path> (or the EFF_SDK_ROOT_DIR environment variable); otherwise ask
# the active Python environment where the package is, which is what a wheel install gives.
if(NOT DEFINED EFF_SDK_ROOT_DIR AND DEFINED ENV{EFF_SDK_ROOT_DIR})
  set(EFF_SDK_ROOT_DIR "$ENV{EFF_SDK_ROOT_DIR}")
endif()
if(NOT DEFINED EFF_SDK_ROOT_DIR)
  execute_process(
    COMMAND python -c "import effcc; print(effcc.__path__[0])"
    OUTPUT_VARIABLE _effcc_pkg
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _effcc_rc)
  if(NOT _effcc_rc EQUAL 0)
    message(
      FATAL_ERROR
        "effcc not found. Activate the Python environment it is installed in, or pass -DEFF_SDK_ROOT_DIR=<effcc package>/sdk"
    )
  endif()
  set(EFF_SDK_ROOT_DIR "${_effcc_pkg}/sdk")
  # Tell the SDK where the compiler binaries are as well. On Windows the environment's
  # effcc launcher lives in Scripts\, which the SDK's PATH lookup does not understand.
  if(NOT DEFINED EFFCC_DIR)
    set(EFFCC_DIR "${_effcc_pkg}")
  endif()
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
