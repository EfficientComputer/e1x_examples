add_eff_app(recursive_factorial
    TYPE exe
    ARCHS e1x
    TARGETS fabric
    SOURCE main.c recursive_factorial.c
)

eff_subtarget_compile_options(recursive_factorial_fabric PRIVATE --tailcall --skip-iterative-mapping)
