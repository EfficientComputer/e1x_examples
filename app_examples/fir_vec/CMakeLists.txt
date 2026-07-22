add_eff_app(fir_vec
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c fir_vec.c
)

eff_subtarget_compile_options(fir_vec_fabric PRIVATE -DUSE_PRAGMAS)
