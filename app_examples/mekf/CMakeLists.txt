add_eff_app(mekf
    TYPE exe
    ARCHS e1x
    TARGETS scalar
    SOURCE main.c mekf_q20.c
)

# mekf_q20.c uses libm (sqrtf in init, acosf in the readout, fmaf from fused
# float multiply-adds) — link it explicitly so it resolves on toolchains that
# don't auto-link libm. (The MEKF hot path is pure int32; float is init-only.)
target_link_libraries(mekf_scalar m)
