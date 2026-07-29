add_eff_app(dmv_vec
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c dmv_vec.c
)

eff_subtarget_compile_options(dmv_vec_fabric PRIVATE -I ${CMAKE_SOURCE_DIR}/shared -DUSE_PRAGMAS)
