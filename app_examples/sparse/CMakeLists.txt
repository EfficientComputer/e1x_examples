add_eff_app(sparse
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c svec.c
)

eff_subtarget_compile_options(sparse_fabric PRIVATE -DUSE_PRAGMAS)
