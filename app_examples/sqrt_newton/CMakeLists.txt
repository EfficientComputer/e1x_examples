add_eff_app(sqrt_newton
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE sqrt_newton.c
)

eff_subtarget_compile_options(sqrt_newton_fabric PRIVATE -DUSE_PRAGMAS)
