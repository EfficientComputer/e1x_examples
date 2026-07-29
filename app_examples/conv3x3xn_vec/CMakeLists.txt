add_eff_app(conv3x3xn_vec
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c conv3x3xn_vec.c
)

eff_subtarget_compile_options(conv3x3xn_vec_fabric PRIVATE -DUSE_PRAGMAS)
