/*
 * F450 sensor layer — fixed-point (see sensors.h). Reads BMI323 + ENS220 and
 * converts to Q-format with pure integer math (no float / no soft-float).
 */
#include "sensors.h"
#include "bmi323.h"
#include "ens220.h"

/* gyro: BMI323 ±2000 °/s -> rad/s in Q15.
 *   rad/s per LSB = (2000/32768)*(pi/180); ×2^15 = 2000*pi/180 = 34.9066.
 *   gyro_q15 = raw * 34907 / 1000  (pure integer, <0.001% error). */
#define GYR_NUM  34907
#define GYR_DEN  1000

/* accel: BMI323 ±2 g (16384 counts/g). mekf wants a_meas/g in Q8:
 *   raw/16384 * 256 = raw/64 = raw >> 6. */
#define ACC_Q8_SHIFT  6

/* baro: ENS220 pressure raw is 1/64 Pa; altitude ~ (P0-P)/(rho*g), rho*g≈12.02.
 *   alt_q20 = (P0_raw - raw_p) * 2^20 / (64 * 12.02) = (P0_raw - raw_p) * 1363. */
#define BARO_ALT_MUL  1363

static int32_t g_p0_raw = 0;

int sensors_init(void)
{
    int mask = 0;
    if (bmi323_init() == 0 && bmi323_configure() == 0) mask |= 0x1;
    if (ens220_init() == 0) {
        uint32_t rp; uint16_t rt;
        if (ens220_read_raw(&rp, &rt) == 0) { g_p0_raw = (int32_t)rp; mask |= 0x2; }
    }
    return mask;
}

void sensors_read(sensors_t *s)
{
    int16_t ax, ay, az, gx, gy, gz;

    if (bmi323_read_gyro(&gx, &gy, &gz) == 0) {
        s->gyro_q15[0] = (int32_t)gx * GYR_NUM / GYR_DEN;
        s->gyro_q15[1] = (int32_t)gy * GYR_NUM / GYR_DEN;
        s->gyro_q15[2] = (int32_t)gz * GYR_NUM / GYR_DEN;
    }
    if (bmi323_read_accel(&ax, &ay, &az) == 0) {
        s->a_norm_q8[0] = (int32_t)ax >> ACC_Q8_SHIFT;
        s->a_norm_q8[1] = (int32_t)ay >> ACC_Q8_SHIFT;
        s->a_norm_q8[2] = (int32_t)az >> ACC_Q8_SHIFT;
    }

    uint32_t rp; uint16_t rt;
    if (ens220_read_raw(&rp, &rt) == 0)
        s->baro_alt_q20 = (g_p0_raw - (int32_t)rp) * BARO_ALT_MUL;
}
