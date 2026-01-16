add_eff_app(ldpc
    TYPE exe
    ARCHS e1x
    TARGETS fabric scalar
    SOURCE generator_80211.c parity_80211.c ldpc.c main.c)

eff_subtarget_compile_options(ldpc_fabric PRIVATE -DUSE_PRAGMAS)

eff_subtarget_link_options(ldpc_scalar PRIVATE -lm)
eff_subtarget_link_options(ldpc_fabric PRIVATE -lm)