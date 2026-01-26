add_eff_app(bme280
    TYPE exe
    ARCHS e1x
    TARGETS scalar
    SOURCE main.c bme280.c
    )

# Enable 32-bit integer mode
# Note: 32-bit mode pressure may need calibration adjustment
target_compile_definitions(bme280_scalar PRIVATE BME280_32BIT_ENABLE)