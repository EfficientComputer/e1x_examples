# @OPT_LEVEL=0
# @OPT_TAGS=
# @OPT_COMMENT=
# @TYPE_TAGS=ml,benchmark
# @OP_COUNTS=1:{529416}
# @AUTOBENCH=all
# @AUTOBENCH_SIZES=1

set(MODEL_NAME ad)

get_filename_component(RIPTOOLS_BINARY_PATH "${CMAKE_C_COMPILER}" DIRECTORY)
set(EFF_IMPORT_PATH ${RIPTOOLS_BINARY_PATH}/eff-import)
set(MODEL_INPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/ad01_int8.tflite)

# --- eff-import: TFLite -> MLIR ---

set(MODEL_OUTPUT_FABRIC_MLIR ${CMAKE_CURRENT_BINARY_DIR}/${MODEL_NAME}.tosa.fabric.mlir)
add_custom_command(
  OUTPUT  ${MODEL_OUTPUT_FABRIC_MLIR}
  COMMAND ${EFF_IMPORT_PATH} ${MODEL_INPUT_PATH}
          -o ${MODEL_OUTPUT_FABRIC_MLIR}
          -S
  DEPENDS ${MODEL_INPUT_PATH}
  VERBATIM
)
add_custom_target(${MODEL_NAME}-tosa-fabric DEPENDS ${MODEL_OUTPUT_FABRIC_MLIR})

set(MODEL_OUTPUT_SCALAR_MLIR ${CMAKE_CURRENT_BINARY_DIR}/${MODEL_NAME}.tosa.scalar.mlir)
add_custom_command(
  OUTPUT  ${MODEL_OUTPUT_SCALAR_MLIR}
  COMMAND ${EFF_IMPORT_PATH} ${MODEL_INPUT_PATH}
          -o ${MODEL_OUTPUT_SCALAR_MLIR}
          --no-rip-region
          -S
  DEPENDS ${MODEL_INPUT_PATH}
  VERBATIM
)
add_custom_target(${MODEL_NAME}-tosa-scalar DEPENDS ${MODEL_OUTPUT_SCALAR_MLIR})

# --- Build ---

add_eff_app(mlperftiny_ad
  TYPE exe
  ARCHS e1x e1
  TARGETS fabric scalar # sim
  SOURCE
    ${CMAKE_CURRENT_SOURCE_DIR}/eff/benchmark.c
    ${CMAKE_CURRENT_SOURCE_DIR}/main.c
)

eff_subtarget_sources(mlperftiny_${MODEL_NAME}_fabric ${MODEL_OUTPUT_FABRIC_MLIR})
eff_subtarget_sources(mlperftiny_${MODEL_NAME}_sim ${MODEL_OUTPUT_FABRIC_MLIR})
set_source_files_properties(${MODEL_OUTPUT_FABRIC_MLIR} PROPERTIES LANGUAGE C)

eff_subtarget_sources(mlperftiny_${MODEL_NAME}_scalar ${MODEL_OUTPUT_SCALAR_MLIR})
set_source_files_properties(${MODEL_OUTPUT_SCALAR_MLIR} PROPERTIES LANGUAGE C)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
if(TARGET mlperftiny_${MODEL_NAME}_fabric)
    add_dependencies(mlperftiny_${MODEL_NAME}_fabric ${MODEL_NAME}-tosa-fabric)
endif()
if(TARGET mlperftiny_${MODEL_NAME}_sim)
    add_dependencies(mlperftiny_${MODEL_NAME}_sim ${MODEL_NAME}-tosa-fabric)
endif()
if(TARGET mlperftiny_${MODEL_NAME}_scalar)
    add_dependencies(mlperftiny_${MODEL_NAME}_scalar ${MODEL_NAME}-tosa-scalar)
endif()
