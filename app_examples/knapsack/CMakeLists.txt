add_eff_app(knapsack
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c knapsack.c
)

eff_subtarget_compile_options(knapsack_fabric PRIVATE -DUSE_PRAGMAS --disable-memory-ordering --mof-passes=convert-parallel-for,lift-scf-parallel --disable-mof-option-exploration)
