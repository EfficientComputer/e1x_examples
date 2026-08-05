#pragma once
#include <stdint.h>

/*
 * F450 sensor layer — FIXED-POINT (no float / no soft-float).
 *
 * Inner-loop subset: BMI323 (gyro+accel) and ENS220 (baro). Outputs are in the
 * exact fixed formats the fixed-point estimator consumes:
 *   gyro_q15   : body rate, rad/s in Q15
 *   a_norm_q8  : specific force / g, body frame, Q8 (=256 at 1 g) — mekf_q20 input
 *   baro_alt_q20: altitude relative to startup, metres in Q20 (+up)
 *
 * (BMM350 mag and NEO-M9N GPS join in the position/yaw fixed-point follow-up.)
 *
 * BRING-UP: verify axis signs/mounting vs the airframe body frame; the BMI323
 * driver configures ±2 g (saturates above hover — use ±8 g for real flight).
 */
typedef struct {
    int32_t gyro_q15[3];     /* rad/s, Q15                     */
    int32_t a_norm_q8[3];    /* specific force / g, Q8         */
    int32_t baro_alt_q20;    /* altitude (m, +up), Q20         */
} sensors_t;

/* Init IMU + baro; capture the baro reference pressure. Returns bit0=imu, bit1=baro. */
int  sensors_init(void);

/* Read one fixed-point sample set. */
void sensors_read(sensors_t *s);
