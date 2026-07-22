add_eff_app(fft
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c fft.c fft.h ${CMAKE_CURRENT_BINARY_DIR}/expected.c ${CMAKE_CURRENT_BINARY_DIR}/fft_twiddles.c
)

eff_subtarget_compile_options(fft_fabric PRIVATE -DUSE_PRAGMAS --disable-memory-ordering)

add_custom_command(
    OUTPUT fft_twiddles.c
           expected.c
    COMMAND python3 generate_fft.py --size 4096 -o ${CMAKE_CURRENT_BINARY_DIR}/fft_twiddles.c
    COMMAND python3 generate_fft.py --size 4096 --expected -o ${CMAKE_CURRENT_BINARY_DIR}/expected.c
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/generate_fft.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
