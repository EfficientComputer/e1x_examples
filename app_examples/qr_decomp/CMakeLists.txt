add_eff_app(qr_decomp
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c matrix.c qr.c qr_decomp.c
)

eff_subtarget_compile_options(qr_decomp_fabric PRIVATE -DUSE_PRAGMAS)
