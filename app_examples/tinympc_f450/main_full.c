/*
 * End-to-end TinyMPC flight stack for a Hawk's Work F450 — FIXED-POINT inner loop.
 *
 *   BMI323 (gyro+accel) ─┐ sensors.c (int)   ┌─ mekf_q20 (attitude, Q15/Q20 int)
 *   ENS220 (baro) ───────┘ ───────────────── ┤  + gyro-integrated yaw + baro z
 *                                             └─> 12-state estimate x̂ (Q20)
 *                                                     │
 *                              e = x̂ − xref(t) ─> TinyMPC (Q20, fabric) ─> z[0..3]
 *                                                     └─> f450_motors -> PCA9685 -> ESCs
 *
 * Everything on the per-step path is int32 fixed point — NO float, NO soft-float.
 * The attitude filter is the validated mekf_q20 (gyro+accel; yaw from gyro
 * integration). Horizontal position/velocity and mag-yaw fusion are the staged
 * fixed-point follow-up, so px,py,vx,vy read 0 here and the geofence is inactive:
 * THIS BUILD IS BENCH / TETHERED ONLY (props off) until the position filter lands.
 *
 * Runs the in-place cage routine (hover / vertical bounce / yaw spin) at 50 Hz.
 */
#include <eff.h>
#include <eff/mtimer.h>
#include <stdio.h>
#include <stdint.h>
#include "tinympc.h"
#include "f450_cache_q20.h"
#include "sensors.h"
#include "mekf_q20.h"
#include "f450_motors.h"

#define ONE       (1 << F450Q)          /* Q20 */
#define HZ        50
#define N_RUN     1500                  /* 30 s bring-up routine */
#define UPDATE_DIV 5                    /* accel update every N predicts */

/* per-step estimator params, precomputed to integer literals (no init float):
 *   dt=0.02 -> dt_q15=655; qg=qbg=1 (min-clamped); Ra≈2 */
#define DT_Q15    655
#define QG_Q20    1
#define QBG_Q20   1
#define RA_Q20    2

static int32_t z[NHORIZON*NU], y[NHORIZON*NU], pp[(NHORIZON+1)*NX], dd[NHORIZON*NU], xx[(NHORIZON+1)*NX], uu[NHORIZON*NU];

/* reference: hover -> vertical bounce -> yaw spin -> hover (Q20, 12-state) */
static void f450_ref(int s, int32_t xr[12])
{
    for (int i = 0; i < 12; i++) xr[i] = 0;
    const int32_t ZAMP = (int32_t)(0.8 * ONE);
    if (s >= 200 && s < 800) {                          /* two triangle bounces */
        int local = (s - 200) % 300;
        if (local < 150) { xr[2] = (int32_t)((int64_t)ZAMP*local/150);       xr[8] =  ZAMP/3; }
        else             { xr[2] = (int32_t)((int64_t)ZAMP*(300-local)/150); xr[8] = -ZAMP/3; }
    }
    if (s >= 900) {                                     /* yaw ramp 0..2pi then hold */
        int ys = s - 900; if (ys > 500) ys = 500;
        xr[5]  = (int32_t)((int64_t)((int32_t)(6.2831853*ONE)) * ys / 500);
        xr[11] = (ys < 500) ? (int32_t)(6.2831853*ONE) / 10 : 0;
    }
}

int main(void)
{
    sleep_ms(2000);
    printf("[f450-full] fixed-point flight stack: sensors(int) -> mekf_q20 -> MPC(Q20) -> PWM\r\n");

    int mask = sensors_init();
    printf("[f450-full] sensors up: IMU=%d BARO=%d\r\n", !!(mask&1), !!(mask&2));

    int motors_ok = (f450_motors_init() == 0);
    printf(motors_ok ? "[f450-full] PCA9685 armed (ESCs @MIN) -- PROPS OFF; geofence inactive (bench only)\r\n"
                     : "[f450-full] PCA9685 init FAILED -- telemetry only\r\n");

    MEKF_Q20 mekf;
    float q0[4] = {1,0,0,0}, bg0[3] = {0,0,0};
    mekf_q20_init(&mekf, q0, bg0, 0.2f, 0.05f);

    sensors_t sv;
    int32_t yaw_q20 = 0, z_prev = 0, vz_q20 = 0;
    int est = 0;
    eff_mtimer_ticks_t t_next = eff_mtimer_uptime_us();

    for (int k = 0; k < N_RUN; k++) {
        /* ── sense (fixed point) ────────────────────────────────────────── */
        sensors_read(&sv);

        /* ── attitude: mekf_q20 (int32 hot path) ───────────────────────────*/
        mekf_q20_predict(&mekf, sv.gyro_q15, DT_Q15, QG_Q20, QBG_Q20);
        if (++est >= UPDATE_DIV) { est = 0; mekf_q20_update_accel(&mekf, sv.a_norm_q8, RA_Q20); }

        /* yaw by integrating bias-corrected gyro-z (Q15 rad/s -> Q20 rad) */
        int32_t wz_c = sv.gyro_q15[2] - mekf.b_g[2];
        yaw_q20 += (int32_t)(((int64_t)wz_c * DT_Q15) >> 15) << 5;

        /* baro altitude + low-passed vertical velocity */
        int32_t zq = sv.baro_alt_q20;
        int32_t vz_raw = (zq - z_prev) * HZ;
        vz_q20 += (vz_raw - vz_q20) >> 2;
        z_prev = zq;

        /* ── assemble 12-state estimate (Q20) ─────────────────────────────
         * px,py,vx,vy = 0 (staged position filter). roll/pitch small-angle
         * from the Q15 quaternion vector part (valid near hover). */
        int32_t x12[12];
        x12[0] = 0;                             x12[1] = 0;
        x12[2] = zq;
        x12[3] = (2 * mekf.q[1]) << 5;          /* roll  */
        x12[4] = (2 * mekf.q[2]) << 5;          /* pitch */
        x12[5] = yaw_q20;                       /* yaw   */
        x12[6] = 0;                             x12[7] = 0;
        x12[8] = vz_q20;
        x12[9]  = (sv.gyro_q15[0] - mekf.b_g[0]) << 5;   /* wx */
        x12[10] = (sv.gyro_q15[1] - mekf.b_g[1]) << 5;   /* wy */
        x12[11] = wz_c << 5;                              /* wz */

        /* ── control: e = x̂ − xref (Q20) ──────────────────────────────────*/
        int32_t xr[12], e[12];
        f450_ref(k, xr);
        for (int i = 0; i < 12; i++) e[i] = x12[i] - xr[i];
        tinympc_solve_fabric(e, F450_A_Q20, F450_B_Q20, F450_Kinf_Q20, F450_Quu_Q20, F450_AmBKt_Q20,
            F450_UMIN_Q20, F450_UMAX_Q20, RHO, NHORIZON, NHORIZON, z, y, pp, dd, xx, uu);

        if (motors_ok) f450_motors_write(z);

        if (k % 100 == 0)
            printf("[f450-full] %4d  z=%ldmm zref=%ldmm  rpy=[%ld %ld %ld]deg  u0=%ld\r\n", k,
                   (long)x12[2]*1000/ONE, (long)xr[2]*1000/ONE,
                   (long)((int64_t)x12[3]*57296/ONE/1000),
                   (long)((int64_t)x12[4]*57296/ONE/1000),
                   (long)((int64_t)x12[5]*57296/ONE/1000),
                   (long)z[0]*1000/ONE);

        t_next += 1000000 / HZ;
        while (eff_mtimer_uptime_us() < t_next) { }
    }

    if (motors_ok) f450_motors_disarm();
    printf("[f450-full] routine complete -- motors disarmed\r\n");
    return 0;
}
