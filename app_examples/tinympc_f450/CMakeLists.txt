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

# ── End-to-end flight stack: real sensors -> fixed-point estimator -> MPC -> PWM ──
# Inner loop is all int32 fixed point: BMI323+ENS220 -> mekf_q20 attitude -> TinyMPC
# Q20 -> PCA9685. (BMM350 mag + NEO-M9N GPS + position fusion are the staged
# fixed-point follow-up.)
add_eff_app(tinympc_f450_full
    TYPE exe
    ARCHS e1x
    TARGETS fabric scalar
    SOURCE main_full.c tinympc.c pca9685.c f450_motors.c sensors.c mekf_q20.c posyaw.c
           bmi323.c bmm350.c ens220.c neo_m9n.c
)

eff_subtarget_compile_options(tinympc_f450_full_fabric PRIVATE -DUSE_PRAGMAS --disable-memory-ordering)
target_compile_definitions(tinympc_f450_full_fabric PRIVATE F450_MOTORS)
target_compile_definitions(tinympc_f450_full_scalar PRIVATE F450_MOTORS)

# The estimator + sensor/quaternion math use libm (sqrtf, atan2f, asinf, cosf);
# link it explicitly so it resolves on toolchains that don't auto-link libm.
target_link_libraries(tinympc_f450_full_fabric m)
target_link_libraries(tinympc_f450_full_scalar m)

# ── Sensor-check / calibration mode: estimator only, NO MOTORS, CSV telemetry ──
add_eff_app(tinympc_f450_check
    TYPE exe
    ARCHS e1x
    TARGETS scalar
    SOURCE main_check.c sensors.c mekf_q20.c posyaw.c
           bmi323.c bmm350.c ens220.c neo_m9n.c
)
target_link_libraries(tinympc_f450_check_scalar m)
