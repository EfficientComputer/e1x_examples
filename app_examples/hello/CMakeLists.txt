add_eff_app(hello
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE hello.c
)

eff_subtarget_compile_options(hello_fabric PRIVATE -DUSE_PRAGMAS)
