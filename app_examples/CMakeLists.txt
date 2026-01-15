cmake_minimum_required(VERSION 3.15)

project(e1x_examples)

# EVK Requires UART peripheral 3
set(EFF_STDIO_PORT 3)

# Replace with your EFF SDK path if installed elsewhere
file(REAL_PATH "$ENV{HOME}/effcc/sdk" EFF_SDK_ROOT_DIR)
include(${EFF_SDK_ROOT_DIR}/setup_sdk.cmake)
add_subdirectory(${EFF_SDK_ROOT_DIR} sdk)

# A NOTE TO PRODUCT TEAM: REMOVE THE FOLLOWING LINE AND
# VALIDATE WITH FINAL RELEASE VERSION OF COMPILER
add_compile_definitions(__efficient__=__effcc_rip)

add_subdirectory(biquad_filter)
add_subdirectory(conv5x5)
add_subdirectory(fft)
add_subdirectory(jpeg)
add_subdirectory(ldpc)
add_subdirectory(quickstart)