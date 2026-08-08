add_eff_app(sort
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c sort.c
)

eff_subtarget_compile_options(sort_fabric PRIVATE -DUSE_PRAGMAS)
