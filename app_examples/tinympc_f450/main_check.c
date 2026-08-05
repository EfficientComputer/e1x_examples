/*
 * F450 sensor-check / calibration mode — E1x, NO MOTORS.
 *
 * Reads all sensors and runs the full fixed-point estimator (mekf_q20 attitude +
 * posyaw position/velocity + mag-yaw), but never touches the PCA9685 — safe to run
 * with the aircraft in hand. Streams one CSV telemetry line per sample over the
 * STDIO UART for the host monitor (f450_monitor.py) or a plain serial terminal.
 *
 * CSV columns (integers):
 *   F450CHK, step,
 *   gyro_mrads[3], anorm_mg[3], mag_raw[3], baro_mm,
 *   gps_fix, gps_x_mm, gps_y_mm,
 *   roll_mdeg, pitch_mdeg, yaw_mdeg,
 *   pos_x_mm, pos_y_mm, pos_z_mm, vel_x_mms, vel_y_mms, vel_z_mms,
 *   bg_x_mrads, bg_y_mrads, bg_z_mrads
 */
#include <eff.h>
#include <eff/mtimer.h>
#include <stdio.h>
#include <stdint.h>
#include "sensors.h"
#include "mekf_q20.h"
#include "posyaw.h"

#define HZ         50
#define UPDATE_DIV 5
#define DT_Q15     655
#define QG_Q20     1
#define QBG_Q20    1
#define RA_Q20     2
#define KYAW       328
#define PI_Q20     3294199
#define TWO_PI_Q20 6588397
#define PRINT_EVERY 3                 /* ~16 Hz telemetry */

static int32_t wrap_pi(int32_t a){ while(a> PI_Q20)a-=TWO_PI_Q20; while(a<-PI_Q20)a+=TWO_PI_Q20; return a; }

/* Q15 rad -> millideg;  Q15 rad/s -> mrad/s;  Q8 (/g) -> milli-g */
static long q15_mdeg(int32_t q){ return (long)((int64_t)q*57296/32768); }
static long q15_mrads(int32_t q){ return (long)((int64_t)q*1000/32768); }
static long q20_mm(int32_t q){ return (long)((int64_t)q*1000/1048576); }
static long q20_mdeg(int32_t q){ return (long)((int64_t)q*57296/1048576); }

int main(void)
{
    sleep_ms(2000);
    int mask = sensors_init();
    printf("# F450 sensor-check (NO MOTORS). sensors: IMU=%d BARO=%d MAG=%d GPS=%d\r\n",
           !!(mask&1), !!(mask&2), !!(mask&4), !!(mask&8));
    printf("# CSV: F450CHK,step,gyro_mrads xyz,anorm_mg xyz,mag xyz,baro_mm,"
           "gps_fix,gpsx_mm,gpsy_mm,roll_md,pitch_md,yaw_md,px_mm,py_mm,pz_mm,vx,vy,vz,bgx,bgy,bgz\r\n");

    MEKF_Q20 mekf; float q0[4]={1,0,0,0}, bg0[3]={0,0,0};
    mekf_q20_init(&mekf, q0, bg0, 0.2f, 0.05f);
    poskf_t pk; poskf_init(&pk);

    sensors_t sv; int32_t yaw_q20 = 0; int est = 0;
    eff_mtimer_ticks_t t_next = eff_mtimer_uptime_us();

    for (int k = 0; ; k++) {
        sensors_read(&sv);
        mekf_q20_predict(&mekf, sv.gyro_q15, DT_Q15, QG_Q20, QBG_Q20);
        if (++est >= UPDATE_DIV) { est = 0; mekf_q20_update_accel(&mekf, sv.a_norm_q8, RA_Q20); }
        poskf_predict(&pk, mekf.q, sv.accel_q20);
        poskf_update_baro(&pk, sv.baro_alt_q20);
        if (sv.gps_valid) poskf_update_gps(&pk, sv.gps_x_q20, sv.gps_y_q20);
        int32_t wz_c = sv.gyro_q15[2] - mekf.b_g[2];
        yaw_q20 += (int32_t)(((int64_t)wz_c * DT_Q15) >> 15) << 5;
        int32_t mag_h = (int32_t)posyaw_mag_heading_q15(mekf.q, sv.mag) << 5;
        yaw_q20 += (int32_t)(((int64_t)KYAW * wrap_pi(mag_h - yaw_q20)) >> 15);

        if (k % PRINT_EVERY == 0) {
            printf("F450CHK,%d,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%d,%ld,%ld,"
                   "%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld\r\n", k,
                q15_mrads(sv.gyro_q15[0]), q15_mrads(sv.gyro_q15[1]), q15_mrads(sv.gyro_q15[2]),
                (long)sv.a_norm_q8[0]*1000/256, (long)sv.a_norm_q8[1]*1000/256, (long)sv.a_norm_q8[2]*1000/256,
                (long)sv.mag[0], (long)sv.mag[1], (long)sv.mag[2],
                q20_mm(sv.baro_alt_q20),
                sv.gps_valid, q20_mm(sv.gps_x_q20), q20_mm(sv.gps_y_q20),
                q20_mdeg((2*mekf.q[1])<<5), q20_mdeg((2*mekf.q[2])<<5), q20_mdeg(yaw_q20),
                q20_mm(pk.r[0]), q20_mm(pk.r[1]), q20_mm(pk.r[2]),
                q20_mm(pk.v[0]), q20_mm(pk.v[1]), q20_mm(pk.v[2]),
                q15_mrads(mekf.b_g[0]), q15_mrads(mekf.b_g[1]), q15_mrads(mekf.b_g[2]));
        }

        t_next += 1000000 / HZ;
        while (eff_mtimer_uptime_us() < t_next) { }
    }
    return 0;
}
