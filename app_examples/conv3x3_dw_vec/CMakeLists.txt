add_eff_app(conv3x3_dw_vec
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c conv3x3_dw_vec.c
)

eff_subtarget_compile_options(conv3x3_dw_vec_fabric PRIVATE -DUSE_PRAGMAS)
