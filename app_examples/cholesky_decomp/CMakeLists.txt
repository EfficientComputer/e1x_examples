add_eff_app(cholesky_decomp
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c cholesky.c
)

eff_subtarget_compile_options(cholesky_decomp_fabric PRIVATE -DUSE_PRAGMAS)

if(EFF_BLD_HAND_OPTIMIZED)
    eff_subtarget_compile_options(cholesky_decomp_fabric PRIVATE --disable-memory-ordering)
endif()
