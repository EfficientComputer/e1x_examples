add_eff_app(conv5x5sym
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c conv5x5sym.c
)

eff_subtarget_compile_options(conv5x5sym_fabric PRIVATE -DUSE_PRAGMAS)
