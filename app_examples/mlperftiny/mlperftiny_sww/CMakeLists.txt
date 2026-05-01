# @OPT_LEVEL=0
# @OPT_TAGS=
# @OPT_COMMENT=
# @TYPE_TAGS=ml,benchmark
# @OP_COUNTS=1:{1703877}
# @AUTOBENCH=all
# @AUTOBENCH_SIZES=1
# @AUTOBENCH_OPTIMIZED_FLAGS=


set(MODEL_NAME sww)

get_filename_component(RIPTOOLS_BINARY_PATH "${CMAKE_C_COMPILER}" DIRECTORY)
set(EFF_IMPORT_PATH ${RIPTOOLS_BINARY_PATH}/eff-import)
set(MODEL_INPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/str_ww_ref_model.tflite)

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
set_source_files_properties(${MODEL_OUTPUT_FABRIC_MLIR} PROPERTIES LANGUAGE C)
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
set_source_files_properties(${MODEL_OUTPUT_SCALAR_MLIR} PROPERTIES LANGUAGE C)
add_custom_target(${MODEL_NAME}-tosa-scalar DEPENDS ${MODEL_OUTPUT_SCALAR_MLIR})

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
set_source_files_properties(${MODEL_OUTPUT_OPT_FABRIC_MLIR} PROPERTIES LANGUAGE C)
add_custom_target(${MODEL_NAME}-tosa-opt-fabric DEPENDS ${MODEL_OUTPUT_OPT_FABRIC_MLIR})

# --- Hand-written MLIR ---

set(MLIR_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/eff/sdot_vec4xi8.mlir
)

# -----------------------------------------------------------------------------
# Final step
# -----------------------------------------------------------------------------

  set(COMMON_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/eff/benchmark.c
    ${CMAKE_CURRENT_SOURCE_DIR}/main.c
    ${CMAKE_CURRENT_SOURCE_DIR}/featurize.c
    ${CMAKE_CURRENT_SOURCE_DIR}/fft4k.c
    ${CMAKE_CURRENT_SOURCE_DIR}/eff/optimized_kernels_${EFF_ARCH}.c
  )

add_eff_app(mlperftiny_sww
  TYPE exe
  ARCHS e1x e1
  TARGETS fabric optimized_fabric scalar # optimized_sim
  SOURCE
    ${COMMON_SRCS}
    hamming_coeffs.h
    fft4_1024.h
    mel_filterbank.h
    ${MLIR_SRCS}
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
if(TARGET mlperftiny_${MODEL_NAME}_fabric)
  eff_subtarget_sources(mlperftiny_${MODEL_NAME}_fabric ${MODEL_OUTPUT_FABRIC_MLIR})
  add_dependencies(mlperftiny_${MODEL_NAME}_fabric      ${MODEL_NAME}-tosa-fabric)
endif()
if(TARGET mlperftiny_${MODEL_NAME}_sim)
  eff_subtarget_sources(mlperftiny_${MODEL_NAME}_sim    ${MODEL_OUTPUT_FABRIC_MLIR})
  add_dependencies(mlperftiny_${MODEL_NAME}_sim         ${MODEL_NAME}-tosa-fabric)
endif()
if(TARGET mlperftiny_${MODEL_NAME}_scalar)
  eff_subtarget_sources(mlperftiny_${MODEL_NAME}_scalar ${MODEL_OUTPUT_SCALAR_MLIR})
  add_dependencies(mlperftiny_${MODEL_NAME}_scalar      ${MODEL_NAME}-tosa-scalar)
endif()
if(TARGET mlperftiny_${MODEL_NAME}_optimized_fabric)
  eff_subtarget_sources(mlperftiny_${MODEL_NAME}_optimized_fabric ${MODEL_OUTPUT_OPT_FABRIC_MLIR})
  add_dependencies(mlperftiny_${MODEL_NAME}_optimized_fabric      ${MODEL_NAME}-tosa-opt-fabric)
  eff_subtarget_compile_options(mlperftiny_${MODEL_NAME}_optimized_fabric PRIVATE -DUSE_PRAGMAS)
endif()
if(TARGET mlperftiny_${MODEL_NAME}_optimized_sim)
  eff_subtarget_sources(mlperftiny_${MODEL_NAME}_optimized_sim    ${MODEL_OUTPUT_OPT_FABRIC_MLIR})
  add_dependencies(mlperftiny_${MODEL_NAME}_optimized_sim         ${MODEL_NAME}-tosa-opt-fabric)
  eff_subtarget_compile_options(mlperftiny_${MODEL_NAME}_optimized_sim    PRIVATE -DUSE_PRAGMAS)
endif()


# Note: the following custom commands are used to generate header files for the application.
# The files are checked into the repository for simplicity, but can be regenerated using these commands if needed.
# add_custom_command(
#   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/hamming_coeffs.h
#   COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && python3 gen_hamming_window.py --window_size 1024 --output ${CMAKE_CURRENT_SOURCE_DIR}/hamming_coeffs.h
# )

# add_custom_command(
#   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/fft4_1024.h
#   COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && python3 generate_fft4.py --size 1024 --output ${CMAKE_CURRENT_SOURCE_DIR}/fft4_1024.h --chip=${EFF_ARCH}
# )

# add_custom_command(
#   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/mel_filterbank.h
#   COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && python3 gen_mel_filterbank.py -o ${CMAKE_CURRENT_SOURCE_DIR}/mel_filterbank.h
# )
