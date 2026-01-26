include_directories(biquad_filter ${CMAKE_CURRENT_BINARY_DIR})

add_eff_app(biquad_filter
    TYPE exe
    ARCHS e1x
    TARGETS fabric scalar
    SOURCE main.c biquad_filter.c biquad_filter_input.h.inc
    )

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/biquad_filter_input.h.inc
  COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && python3 gen_data.py -n 16384 -o ${CMAKE_CURRENT_BINARY_DIR}/biquad_filter_input.h.inc --freq1=200 --freq2=2000 --low-gain=0.86 --phase-shift=-1
)