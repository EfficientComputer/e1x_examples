add_eff_app(crc32
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c crc32.c
)

eff_subtarget_compile_options(crc32_fabric PRIVATE -DUSE_PRAGMAS)
