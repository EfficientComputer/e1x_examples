/*
 * F450 sensor layer — fixed-point (see sensors.h). Integer scaling only; the sole
 * float is a one-time cosf when the GPS local-tangent origin is captured.
 */
#include "sensors.h"
#include "bmi323.h"
#include "bmm350.h"
#include "ens220.h"
#include "neo_m9n.h"
#include <math.h>

#define GYR_NUM  34907        /* gyro raw -> Q15 rad/s: raw*34907/1000            */
#define GYR_DEN  1000
#define ACC_Q8_SHIFT 6        /* accel raw -> a_meas/g Q8: raw>>6 (±2 g)          */
#define ACC_Q20_MUL  628      /* accel raw -> m/s² Q20: raw*628                    */
#define BARO_ALT_MUL 1363     /* (P0-P)raw -> altitude Q20                          */
#define GPS_NORTH_MUL 11660   /* (lat_1e7 - lat0) -> north metres Q20              */
#define GPS_RATE_HZ  5

static int32_t g_p0_raw = 0;
static int     g_gps_org = 0;
static int32_t g_lat0, g_lon0, g_east_mul = GPS_NORTH_MUL;

int sensors_init(void)
{
    int mask = 0;
    if (bmi323_init() == 0 && bmi323_configure() == 0) mask |= 0x1;
    if (ens220_init() == 0) {
        uint32_t rp; uint16_t rt;
        if (ens220_read_raw(&rp, &rt) == 0) { g_p0_raw = (int32_t)rp; mask |= 0x2; }
    }
    if (bmm350_init() == 0 && bmm350_configure() == 0) mask |= 0x4;
    if (neo_m9n_init() == 0 && neo_m9n_enable_ubx_pvt(GPS_RATE_HZ) == 0) mask |= 0x8;
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
        s->accel_q20[0] = (int32_t)ax * ACC_Q20_MUL;
        s->accel_q20[1] = (int32_t)ay * ACC_Q20_MUL;
        s->accel_q20[2] = (int32_t)az * ACC_Q20_MUL;
    }

    int32_t mx, my, mz;
    if (bmm350_read_mag(&mx, &my, &mz) == 0) { s->mag[0]=mx; s->mag[1]=my; s->mag[2]=mz; }

    uint32_t rp; uint16_t rt;
    if (ens220_read_raw(&rp, &rt) == 0)
        s->baro_alt_q20 = (g_p0_raw - (int32_t)rp) * BARO_ALT_MUL;

    s->gps_valid = 0;
    neo_m9n_pvt_t pvt;
    if (neo_m9n_read_pvt(&pvt) == 0 && pvt.fix_type >= 3) {
        if (!g_gps_org) {
            g_lat0 = pvt.lat_1e7; g_lon0 = pvt.lon_1e7;
            /* one-time: east scale shrinks by cos(lat0) */
            g_east_mul = (int32_t)(GPS_NORTH_MUL * cosf((float)pvt.lat_1e7 * 1e-7f * 0.0174533f));
            g_gps_org = 1;
        }
        s->gps_x_q20 = (pvt.lat_1e7 - g_lat0) * GPS_NORTH_MUL;   /* north */
        s->gps_y_q20 = (pvt.lon_1e7 - g_lon0) * g_east_mul;      /* east  */
        s->gps_valid = 1;
    }
}
