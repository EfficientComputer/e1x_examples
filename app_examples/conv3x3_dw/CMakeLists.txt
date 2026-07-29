add_eff_app(conv3x3_dw
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c conv3x3_dw.c
)

eff_subtarget_compile_options(conv3x3_dw_fabric PRIVATE -DUSE_PRAGMAS --disable-memory-ordering)
