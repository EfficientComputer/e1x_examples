add_eff_app(fir
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c fir.c
)

eff_subtarget_compile_options(fir_fabric PRIVATE -DUSE_PRAGMAS --disable-memory-ordering --pnr-rounds=10)
