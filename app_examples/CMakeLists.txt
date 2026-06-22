cmake_minimum_required(VERSION 3.15)

project(app_examples)

# EVK Requires UART peripheral 3
set(EFF_STDIO_PORT 3)

# Replace with your EFF SDK path if installed elsewhere
file(REAL_PATH "$ENV{HOME}/effcc/sdk" EFF_SDK_ROOT_DIR)
include(${EFF_SDK_ROOT_DIR}/setup_sdk.cmake)
add_subdirectory(${EFF_SDK_ROOT_DIR} sdk)

add_subdirectory(biquad_filter)
add_subdirectory(conv5x5)
add_subdirectory(dmm)
add_subdirectory(fft)
add_subdirectory(jpeg)
add_subdirectory(ldpc)
add_subdirectory(quickstart)
add_subdirectory(tinympc)
add_subdirectory(tinympc_lpv)