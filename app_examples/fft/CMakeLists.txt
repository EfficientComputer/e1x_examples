add_eff_app(fft
    TYPE exe
    ARCHS e1x
    TARGETS fabric sim scalar
    SOURCE main.c fft.c expected.c
    )

eff_subtarget_compile_options(fft_fabric PRIVATE -DUSE_PRAGMAS --disable-memory-ordering)
eff_subtarget_compile_options(fft_sim PRIVATE -DUSE_PRAGMAS --disable-memory-ordering)