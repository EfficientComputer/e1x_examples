add_eff_app(dijkstras
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c dijkstras.c graph.c
)

eff_subtarget_compile_options(dijkstras_fabric PRIVATE -DUSE_PRAGMAS)
