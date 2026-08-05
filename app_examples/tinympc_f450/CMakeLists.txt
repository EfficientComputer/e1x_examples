add_eff_app(tinympc_f450
    TYPE exe
    ARCHS e1x
    TARGETS fabric scalar
    SOURCE main.c tinympc.c pca9685.c f450_motors.c
)

eff_subtarget_compile_options(tinympc_f450_fabric PRIVATE -DUSE_PRAGMAS --disable-memory-ordering)

# Emit the 4 ESC PWM signals via the PCA9685 (I2C). Drop this to run the
# control loop as a pure sim with no motor output.
target_compile_definitions(tinympc_f450_fabric PRIVATE F450_MOTORS)
target_compile_definitions(tinympc_f450_scalar PRIVATE F450_MOTORS)
