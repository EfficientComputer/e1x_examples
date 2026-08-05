#pragma once
#include <stdint.h>

/*
 * F450 sensor layer — FIXED-POINT (no float in the hot path).
 *
 *   BMI323  gyro+accel   -> gyro_q15, a_norm_q8 (mekf), accel_q20 (position predict)
 *   BMM350  magnetometer -> mag (raw counts; direction only, used via CORDIC)
 *   ENS220  barometer    -> baro_alt_q20
 *   NEO-M9N GPS          -> local north/east (Q20 m) on each 3D fix
 *
 * Frame: body x-fwd/y-right/z-down as mounted (VERIFY signs); world z-up,
 * x-north, y-east. BMI323 is ±2 g here (use ±8 g for real flight).
 * Only one-time float: cosf at GPS-origin capture.
 */
typedef struct {
    int32_t gyro_q15[3];     /* body rate, rad/s Q15                */
    int32_t a_norm_q8[3];    /* specific force / g, Q8 (mekf input) */
    int32_t accel_q20[3];    /* body specific force, m/s² Q20       */
    int32_t mag[3];          /* body field, raw counts (direction)  */
    int32_t baro_alt_q20;    /* altitude rel. startup, m Q20 (+up)  */
    int     gps_valid;       /* 1 = fresh 3D fix this read          */
    int32_t gps_x_q20;       /* local north, m Q20                  */
    int32_t gps_y_q20;       /* local east,  m Q20                  */
} sensors_t;

/* Init all sensors. Returns bitmask: bit0 imu, bit1 baro, bit2 mag, bit3 gps. */
int  sensors_init(void);

/* Read one fixed-point sample set (GPS polled; gps_valid set only on new fix). */
void sensors_read(sensors_t *s);
