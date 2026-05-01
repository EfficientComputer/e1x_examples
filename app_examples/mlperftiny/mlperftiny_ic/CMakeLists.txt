# @OPT_LEVEL=0
# @OPT_TAGS=
# @OPT_COMMENT=
# @TYPE_TAGS=ml,benchmark
# @OP_COUNTS=1:{25113821}
# @AUTOBENCH=all
# @AUTOBENCH_SIZES=1
# @AUTOBENCH_OPTIMIZED_FLAGS=

set(MODEL_NAME ic)

get_filename_component(RIPTOOLS_BINARY_PATH "${CMAKE_C_COMPILER}" DIRECTORY)
set(EFF_IMPORT_PATH ${RIPTOOLS_BINARY_PATH}/eff-import)
set(MODEL_INPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/models/v1.1/pretrainedResnet_quant.tflite)
# set(MODEL_INPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/models/v1.4-pre/pretrainedResnet_quant.tflite)

set(RULES_PATH ${CMAKE_CURRENT_SOURCE_DIR}/eff/tosa_optimized_rules_${EFF_ARCH}.pdl)

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

set(MODEL_OUTPUT_OPT_FABRIC_MLIR ${CMAKE_CURRENT_BINARY_DIR}/${MODEL_NAME}.tosa.opt.fabric.mlir)
add_custom_command(
  OUTPUT  ${MODEL_OUTPUT_OPT_FABRIC_MLIR}
  COMMAND ${EFF_IMPORT_PATH} ${MODEL_INPUT_PATH}
          -o ${MODEL_OUTPUT_OPT_FABRIC_MLIR}
          --use-custom-kernels
          --custom-rule-path=${RULES_PATH}
          --entrypoint-name=run_model_opt
          -S
  DEPENDS ${MODEL_INPUT_PATH} ${RULES_PATH}
  VERBATIM
)
add_custom_target(${MODEL_NAME}-tosa-opt-fabric DEPENDS ${MODEL_OUTPUT_OPT_FABRIC_MLIR})

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

# --- Hand-written MLIR ---

set(MLIR_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/sdot_vec4xi8.mlir
  ${CMAKE_CURRENT_SOURCE_DIR}/add_tree.mlir
  ${CMAKE_CURRENT_SOURCE_DIR}/load_guard.mlir
)

# --- Build ---

set(COMMON_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/eff/benchmark.c
  ${CMAKE_CURRENT_SOURCE_DIR}/main.c
  ${CMAKE_CURRENT_SOURCE_DIR}/eff/optimized_kernels_${EFF_ARCH}.c
)

add_eff_app(mlperftiny_ic
  TYPE exe
  ARCHS e1x e1
  TARGETS fabric optimized_fabric scalar # sim optimized_sim
  SOURCE
    ${COMMON_SRCS}
    ${MLIR_SRCS}
)

eff_subtarget_sources(mlperftiny_${MODEL_NAME}_fabric ${MODEL_OUTPUT_FABRIC_MLIR})
eff_subtarget_sources(mlperftiny_${MODEL_NAME}_sim    ${MODEL_OUTPUT_FABRIC_MLIR})
set_source_files_properties(${MODEL_OUTPUT_FABRIC_MLIR} PROPERTIES LANGUAGE C)
eff_subtarget_sources(mlperftiny_${MODEL_NAME}_optimized_fabric ${MODEL_OUTPUT_OPT_FABRIC_MLIR})
eff_subtarget_sources(mlperftiny_${MODEL_NAME}_optimized_sim    ${MODEL_OUTPUT_OPT_FABRIC_MLIR})
set_source_files_properties(${MODEL_OUTPUT_OPT_FABRIC_MLIR} PROPERTIES LANGUAGE C)
eff_subtarget_sources(mlperftiny_${MODEL_NAME}_scalar ${MODEL_OUTPUT_SCALAR_MLIR})
set_source_files_properties(${MODEL_OUTPUT_SCALAR_MLIR} PROPERTIES LANGUAGE C)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
if(TARGET mlperftiny_${MODEL_NAME}_fabric)
  add_dependencies(mlperftiny_${MODEL_NAME}_fabric ${MODEL_NAME}-tosa-fabric)
endif()
if(TARGET mlperftiny_${MODEL_NAME}_sim)
  add_dependencies(mlperftiny_${MODEL_NAME}_sim    ${MODEL_NAME}-tosa-fabric)
endif()
if(TARGET mlperftiny_${MODEL_NAME}_optimized_fabric)
  add_dependencies(mlperftiny_${MODEL_NAME}_optimized_fabric ${MODEL_NAME}-tosa-opt-fabric)
endif()
if(TARGET mlperftiny_${MODEL_NAME}_optimized_sim)
  add_dependencies(mlperftiny_${MODEL_NAME}_optimized_sim    ${MODEL_NAME}-tosa-opt-fabric)
endif()
if(TARGET mlperftiny_${MODEL_NAME}_scalar)
  add_dependencies(mlperftiny_${MODEL_NAME}_scalar ${MODEL_NAME}-tosa-scalar)
endif()

eff_subtarget_compile_options(mlperftiny_${MODEL_NAME}_optimized_fabric PRIVATE -DUSE_PRAGMAS)
eff_subtarget_compile_options(mlperftiny_${MODEL_NAME}_optimized_sim    PRIVATE -DUSE_PRAGMAS)
