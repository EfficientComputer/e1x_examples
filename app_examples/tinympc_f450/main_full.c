/*
 * End-to-end TinyMPC flight stack for a Hawk's Work F450 — FULL FIXED-POINT.
 *
 *   BMI323 + BMM350 + ENS220 + NEO-M9N ─(int)→ estimator (all Q-format):
 *        mekf_q20 attitude · posyaw position/velocity (GPS+baro+accel) ·
 *        gyro-integrated yaw corrected by CORDIC mag heading
 *                              └─> 12-state x̂ (Q20)
 *                                       │
 *          e = x̂ − xref(t) ─→ TinyMPC (Q20, fabric) ─→ f450_motors → PCA9685
 *
 * Every per-step op is int32 fixed point — no float / no soft-float (only float is
 * one-time init: mekf_q20_init and the GPS-origin cosf). Runs the in-place cage
 * routine (hover / vertical bounce / yaw spin) at 50 Hz with a live geofence.
 *
 * SAFETY: arms the ESCs and drives live PWM. Bench-test with PROPS OFF; calibrate
 * HOVER_US/GAIN_US, mag, GPS origin, and the accel range before flight.
 */
#include <eff.h>
#include <eff/mtimer.h>
#include <stdio.h>
#include <stdint.h>
#include "tinympc.h"
#include "f450_cache_q20.h"
#include "sensors.h"
#include "mekf_q20.h"
#include "posyaw.h"
#include "f450_motors.h"

#define ONE        (1 << F450Q)
#define HZ         50
#define N_RUN      1500
#define UPDATE_DIV 5
#define DT_Q15     655
#define QG_Q20     1
#define QBG_Q20    1
#define RA_Q20     2
#define KYAW       328               /* mag-yaw complementary gain, Q15 (~0.01) */
#define PI_Q20     3294199
#define TWO_PI_Q20 6588397
#define GEOFENCE_Q20 1363149         /* 1.3 m */

static int32_t z[NHORIZON*NU], y[NHORIZON*NU], pp[(NHORIZON+1)*NX], dd[NHORIZON*NU], xx[(NHORIZON+1)*NX], uu[NHORIZON*NU];

static int32_t wrap_pi(int32_t a) { while (a >  PI_Q20) a -= TWO_PI_Q20;
                                    while (a < -PI_Q20) a += TWO_PI_Q20; return a; }

/* reference: hover -> vertical bounce -> yaw spin -> hover (Q20, 12-state) */
static void f450_ref(int s, int32_t xr[12])
{
    for (int i = 0; i < 12; i++) xr[i] = 0;
    const int32_t ZAMP = (int32_t)(0.8 * ONE);
    if (s >= 200 && s < 800) {
        int local = (s - 200) % 300;
        if (local < 150) { xr[2] = (int32_t)((int64_t)ZAMP*local/150);       xr[8] =  ZAMP/3; }
        else             { xr[2] = (int32_t)((int64_t)ZAMP*(300-local)/150); xr[8] = -ZAMP/3; }
    }
    if (s >= 900) {
        int ys = s - 900; if (ys > 500) ys = 500;
        xr[5]  = (int32_t)((int64_t)((int32_t)(6.2831853*ONE)) * ys / 500);
        xr[11] = (ys < 500) ? (int32_t)(6.2831853*ONE) / 10 : 0;
    }
}

int main(void)
{
    sleep_ms(2000);
    printf("[f450-full] full fixed-point flight stack (sensors -> mekf_q20+posyaw -> MPC -> PWM)\r\n");

    int mask = sensors_init();
    printf("[f450-full] sensors: IMU=%d BARO=%d MAG=%d GPS=%d\r\n",
           !!(mask&1), !!(mask&2), !!(mask&4), !!(mask&8));

    int motors_ok = (f450_motors_init() == 0);
    printf(motors_ok ? "[f450-full] PCA9685 armed (ESCs @MIN) -- PROPS OFF to bench-test\r\n"
                     : "[f450-full] PCA9685 init FAILED -- telemetry only\r\n");

    MEKF_Q20 mekf;
    float q0[4] = {1,0,0,0}, bg0[3] = {0,0,0};
    mekf_q20_init(&mekf, q0, bg0, 0.2f, 0.05f);
    poskf_t pk; poskf_init(&pk);

    sensors_t sv;
    int32_t yaw_q20 = 0;
    int est = 0, geofence_hit = -1;
    eff_mtimer_ticks_t t_next = eff_mtimer_uptime_us();

    for (int k = 0; k < N_RUN; k++) {
        sensors_read(&sv);

        /* attitude (mekf_q20, int32) */
        mekf_q20_predict(&mekf, sv.gyro_q15, DT_Q15, QG_Q20, QBG_Q20);
        if (++est >= UPDATE_DIV) { est = 0; mekf_q20_update_accel(&mekf, sv.a_norm_q8, RA_Q20); }

        /* position/velocity (posyaw complementary: accel predict + baro/GPS correct) */
        poskf_predict(&pk, mekf.q, sv.accel_q20);
        poskf_update_baro(&pk, sv.baro_alt_q20);
        if (sv.gps_valid) poskf_update_gps(&pk, sv.gps_x_q20, sv.gps_y_q20);

        /* yaw: gyro integrate, corrected toward CORDIC mag heading (drift only) */
        int32_t wz_c = sv.gyro_q15[2] - mekf.b_g[2];
        yaw_q20 += (int32_t)(((int64_t)wz_c * DT_Q15) >> 15) << 5;
        int32_t mag_h = (int32_t)posyaw_mag_heading_q15(mekf.q, sv.mag) << 5;
        yaw_q20 += (int32_t)(((int64_t)KYAW * wrap_pi(mag_h - yaw_q20)) >> 15);

        /* 12-state estimate (Q20) */
        int32_t x12[12];
        x12[0]=pk.r[0]; x12[1]=pk.r[1]; x12[2]=pk.r[2];
        x12[3]=(2*mekf.q[1])<<5; x12[4]=(2*mekf.q[2])<<5; x12[5]=yaw_q20;
        x12[6]=pk.v[0]; x12[7]=pk.v[1]; x12[8]=pk.v[2];
        x12[9]=(sv.gyro_q15[0]-mekf.b_g[0])<<5; x12[10]=(sv.gyro_q15[1]-mekf.b_g[1])<<5; x12[11]=wz_c<<5;

        /* geofence on estimated horizontal position (now live) */
        int32_t axp = x12[0]<0?-x12[0]:x12[0], ayp = x12[1]<0?-x12[1]:x12[1];
        if (axp > GEOFENCE_Q20 || ayp > GEOFENCE_Q20) {
            if (motors_ok) f450_motors_disarm();
            geofence_hit = k;
            printf("[f450-full] GEOFENCE ABORT at step %d\r\n", k);
            break;
        }

        /* control: e = x̂ − xref (Q20) */
        int32_t xr[12], e[12];
        f450_ref(k, xr);
        for (int i = 0; i < 12; i++) e[i] = x12[i] - xr[i];
        tinympc_solve_fabric(e, F450_A_Q20, F450_B_Q20, F450_Kinf_Q20, F450_Quu_Q20, F450_AmBKt_Q20,
            F450_UMIN_Q20, F450_UMAX_Q20, RHO, NHORIZON, NHORIZON, z, y, pp, dd, xx, uu);
        if (motors_ok) f450_motors_write(z);

        if (k % 100 == 0)
            printf("[f450-full] %4d  pos[%ld %ld %ld]mm  yaw=%ldmdeg  u0=%ld\r\n", k,
                   (long)x12[0]*1000/ONE, (long)x12[1]*1000/ONE, (long)x12[2]*1000/ONE,
                   (long)((int64_t)x12[5]*57296/ONE), (long)z[0]*1000/ONE);

        t_next += 1000000 / HZ;
        while (eff_mtimer_uptime_us() < t_next) { }
    }

    if (motors_ok) f450_motors_disarm();
    printf(geofence_hit<0 ? "[f450-full] routine complete -- motors disarmed\r\n"
                          : "[f450-full] aborted -- motors disarmed\r\n");
    return 0;
}
