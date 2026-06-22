add_eff_app(tinympc_lpv
    TYPE exe
    ARCHS e1x
    TARGETS fabric scalar
    SOURCE main.c tinympc_lpv.c
)

eff_subtarget_compile_options(tinympc_lpv_fabric PRIVATE -DUSE_PRAGMAS --disable-memory-ordering)
