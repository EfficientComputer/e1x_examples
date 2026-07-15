add_eff_app(jpeg
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c jpeg.c image.c
)

eff_subtarget_compile_options(jpeg_fabric PRIVATE -DUSE_PRAGMAS --disable-memory-ordering --mof-passes=convert-parallel-for,lift-scf-parallel --disable-mof-option-exploration --mark-critical-paths)
