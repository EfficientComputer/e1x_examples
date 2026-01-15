add_eff_app(dmm
    TYPE exe
    ARCHS e1x
    TARGETS fabric sim scalar
    SOURCE main.c dmm.c
)

eff_subtarget_compile_options(dmm_fabric PRIVATE -DUSE_PRAGMAS --disable-memory-ordering)
eff_subtarget_compile_options(dmm_sim PRIVATE -DUSE_PRAGMAS --disable-memory-ordering)
