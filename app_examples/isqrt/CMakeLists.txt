add_eff_app(isqrt
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c isqrt.c
)

eff_subtarget_compile_options(isqrt_fabric PRIVATE -DUSE_PRAGMAS)
