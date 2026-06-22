add_eff_app(tinympc
    TYPE exe
    ARCHS e1x
    TARGETS fabric scalar
    SOURCE main.c tinympc.c
)

eff_subtarget_compile_options(tinympc_fabric PRIVATE -DUSE_PRAGMAS --disable-memory-ordering)
