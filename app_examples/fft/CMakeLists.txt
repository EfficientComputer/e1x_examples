add_eff_app(fft
    TYPE exe
    ARCHS e1x
    TARGETS fabric scalar
    SOURCE main.c fft.c expected.c
    )

eff_subtarget_compile_options(fft_fabric PRIVATE -DUSE_PRAGMAS --disable-memory-ordering)
