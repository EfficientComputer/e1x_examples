cmake_minimum_required(VERSION 3.15)

project(e1x_sensor_examples)

file(REAL_PATH "~/effcc/sdk" EFF_SDK_ROOT_DIR)
include(${EFF_SDK_ROOT_DIR}/setup_sdk.cmake)
add_subdirectory(${EFF_SDK_ROOT_DIR} sdk)

add_subdirectory(bme280)
add_subdirectory(adxl372)