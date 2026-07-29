add_eff_app(levenshtein
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c levenshtein.c
)

eff_subtarget_compile_options(levenshtein_fabric PRIVATE -DUSE_PRAGMAS)

if(EFF_BLD_HAND_OPTIMIZED)
    eff_subtarget_compile_options(levenshtein_fabric PRIVATE --disable-memory-ordering --pnr-rounds=8 --mof-passes=convert-parallel-for,lift-scf-parallel --disable-mof-option-exploration)
endif()
