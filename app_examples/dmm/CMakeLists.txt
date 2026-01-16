add_eff_app(dmm
    TYPE exe
    ARCHS e1x
    TARGETS fabric scalar
    SOURCE main.c dmm.c
)

eff_subtarget_compile_options(dmm_fabric PRIVATE -DUSE_PRAGMAS --disable-memory-ordering)
