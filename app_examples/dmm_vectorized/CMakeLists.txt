add_eff_app(dmm_vectorized
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c dmm_vectorized.c
)

eff_subtarget_compile_options(dmm_vectorized_fabric PRIVATE -I ${CMAKE_SOURCE_DIR}/shared -DUSE_PRAGMAS)
