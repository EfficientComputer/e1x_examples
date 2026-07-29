add_eff_app(dmv_fp
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c dmv_fp.c
)

eff_subtarget_compile_options(dmv_fp_fabric PRIVATE -DUSE_PRAGMAS --disable-memory-ordering)
