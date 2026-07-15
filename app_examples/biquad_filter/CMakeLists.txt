add_eff_app(biquad_filter
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c biquad_filter.c biquad_filter.h ${CMAKE_CURRENT_BINARY_DIR}/biquad_filter_input.c
)

eff_subtarget_compile_options(biquad_filter_fabric PRIVATE -DUSE_PRAGMAS)

add_custom_command(
    OUTPUT biquad_filter_input.c
    COMMAND python3 gen_data.py -n 16384 --freq1=200 --freq2=2000 --low-gain=0.86 --phase-shift=-1 -o ${CMAKE_CURRENT_BINARY_DIR}/biquad_filter_input.c
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/gen_data.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
