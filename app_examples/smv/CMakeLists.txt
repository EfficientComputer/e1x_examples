add_eff_app(smv
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c smv.c
)

eff_subtarget_compile_options(smv_fabric PRIVATE -I ${CMAKE_SOURCE_DIR}/shared -DUSE_PRAGMAS --disable-memory-ordering --mof-passes=convert-parallel-for,lift-scf-parallel --disable-mof-option-exploration)
