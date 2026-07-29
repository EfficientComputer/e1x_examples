add_eff_app(conv5x5_vec
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c conv5x5_vec.c
)

eff_subtarget_compile_options(conv5x5_vec_fabric PRIVATE -DUSE_PRAGMAS)
