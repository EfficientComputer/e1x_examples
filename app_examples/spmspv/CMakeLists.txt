add_eff_app(spmspv
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c spmspv.c
)

eff_subtarget_compile_options(spmspv_fabric PRIVATE -I ${CMAKE_SOURCE_DIR}/shared -DUSE_PRAGMAS --disable-memory-ordering --pnr-rounds=8 --mof-passes=convert-parallel-for,lift-scf-parallel --disable-mof-option-exploration)
