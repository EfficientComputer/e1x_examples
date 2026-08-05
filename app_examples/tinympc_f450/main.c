/*
 * TinyMPC "in-place" cage demo for a Hawk's Work F450 — E1x demo + check.
 *
 * A wall-safe routine for a tight (~4x4x3 m) cage on IMU+mag+baro+GPS: it only
 * commands motions the sensors resolve well and that never translate toward a
 * wall — a vertical bounce (baro-precise) and a yaw spin-in-place — bracketed by
 * hover holds. Horizontal position is NEVER commanded away from centre; a hard
 * GEOFENCE aborts if the craft drifts past a safe radius (the real-flight guard
 * against GPS wander / disturbance).
 *
 * This is trajectory *tracking*: each 50 Hz step we form the error e = x - xref(t)
 * and regulate it to zero with the TinyMPC ADMM solve on the fabric (Q20). The
 * reference is piecewise-linear (triangle bounce, ramped yaw) so it needs no
 * trig and no float. Cache: f450_cache_q20.h (gen_f450_cache.py).
 */
#include <eff.h>
#include <stdio.h>
#include <stdint.h>
#include "tinympc.h"
#include "f450_cache_q20.h"

#define ONE (1 << F450Q)
#define N_SIM     1500          /* 30 s at 50 Hz */
#define GEOFENCE  ((int32_t)(1.3 * ONE))   /* abort if |px| or |py| exceeds this */
#define YAW_2PI   ((int32_t)(6.2831853 * ONE))

/* scratch (two sets so we can compare fabric vs scalar on the first solve) */
static int32_t z0[NHORIZON*NU], y0[NHORIZON*NU], p0[(NHORIZON+1)*NX], d0[NHORIZON*NU], x0[(NHORIZON+1)*NX], u0[NHORIZON*NU];
static int32_t z1[NHORIZON*NU], y1[NHORIZON*NU], p1[(NHORIZON+1)*NX], d1[NHORIZON*NU], x1[(NHORIZON+1)*NX], u1[NHORIZON*NU];

/* Reference trajectory xref(step), Q20. Only z (bounce) and yaw (spin) move;
 * px,py stay 0 the whole time so nothing is ever commanded toward a wall. */
static void f450_ref(int s, int32_t xr[NX])
{
    for (int i = 0; i < NX; i++) xr[i] = 0;
    const int32_t ZAMP = (int32_t)(0.8 * ONE);        /* +-0.8 m vertical bounce */

    if (s >= 200 && s < 800) {                         /* 12 s: two triangle bounces */
        int local = (s - 200) % 300;                   /* 6 s cycle: 3 s up, 3 s down */
        if (local < 150) { xr[2] = (int32_t)((int64_t)ZAMP*local/150);        xr[8] =  ZAMP/3; }
        else             { xr[2] = (int32_t)((int64_t)ZAMP*(300-local)/150);  xr[8] = -ZAMP/3; }
    }
    if (s >= 900) {                                    /* 10 s: yaw ramp 0..2pi, then hold */
        int ys = s - 900; if (ys > 500) ys = 500;
        xr[5]  = (int32_t)((int64_t)YAW_2PI * ys / 500);
        xr[11] = (ys < 500) ? YAW_2PI / 10 : 0;        /* wz during the ramp (0.628 rad/s) */
    }
}

/* one linear-model step: x+ = A x + B u   (Q20) */
static void model_step(int32_t x[NX], const int32_t u[NU])
{
    int32_t xn[NX];
    for (int i = 0; i < NX; i++) {
        int64_t s = 0;
        for (int j = 0; j < NX; j++) s += (int64_t)F450_A_Q20[i*NX+j] * x[j];
        for (int k = 0; k < NU; k++) s += (int64_t)F450_B_Q20[i*NU+k] * u[k];
        xn[i] = (int32_t)(s >> F450Q);
    }
    for (int i = 0; i < NX; i++) x[i] = xn[i];
}

int main(void)
{
    int32_t x[NX] = {0};        /* start at hover, centred */
    int32_t e[NX];

    printf("[f450-cage] in-place demo: hover -> vertical bounce -> yaw spin -> hover\r\n");
    printf("[f450-cage] geofence +-%ld mm; no horizontal translation commanded\r\n", (long)GEOFENCE*1000/ONE);
    printf("[f450-cage] step  phase        z(mm) zref(mm)  yaw(deg)  horiz|off|(mm)\r\n");

    int ok = 1, geofence_hit = -1;
    long max_ztrack = 0;
    for (int k = 0; k < N_SIM; k++) {
        int32_t xr[NX];
        f450_ref(k, xr);
        for (int i = 0; i < NX; i++) e[i] = x[i] - xr[i];      /* tracking error state */

        tinympc_solve_fabric(e, F450_A_Q20, F450_B_Q20, F450_Kinf_Q20, F450_Quu_Q20, F450_AmBKt_Q20,
            F450_UMIN_Q20, F450_UMAX_Q20, RHO, NHORIZON, NHORIZON, z0, y0, p0, d0, x0, u0);
        if (k == 0) {   /* correctness: fabric must match the scalar reference exactly */
            tinympc_solve_scalar(e, F450_A_Q20, F450_B_Q20, F450_Kinf_Q20, F450_Quu_Q20, F450_AmBKt_Q20,
                F450_UMIN_Q20, F450_UMAX_Q20, RHO, NHORIZON, NHORIZON, z1, y1, p1, d1, x1, u1);
            for (int i = 0; i < NHORIZON*NU; i++) if (z0[i] != z1[i]) ok = 0;
        }

        /* hard geofence: abort if we drift past the safe horizontal radius */
        long ax = (long)(x[0] < 0 ? -x[0] : x[0]);
        long ay = (long)(x[1] < 0 ? -x[1] : x[1]);
        if ((ax > GEOFENCE || ay > GEOFENCE) && geofence_hit < 0) geofence_hit = k;

        long zt = (long)(x[2]-xr[2]); if (zt<0) zt=-zt; zt = zt*1000/ONE; if (zt>max_ztrack) max_ztrack=zt;

        if (k % 150 == 0) {
            const char *ph = (k<200)?"hover":(k<800)?"bounce":(k<900)?"hover":(k<1400)?"yaw-spin":"hover";
            long hoff = (ax>ay?ax:ay)*1000/ONE;
            long yaw_deg = (long)((int64_t)x[5] * 57296 / ONE / 1000);   /* rad->deg, int64-safe */
            printf("[f450-cage] %4d  %-9s  %6ld  %6ld    %7ld       %6ld\r\n", k, ph,
                   (long)x[2]*1000/ONE, (long)xr[2]*1000/ONE, yaw_deg, hoff);
        }
        if (geofence_hit >= 0) break;
        model_step(x, z0);
    }

    if (geofence_hit >= 0) {
        printf("[f450-cage] GEOFENCE ABORT at step %d (drifted past safe radius)\r\n", geofence_hit);
        ok = 0;
    } else {
        long fz=(long)x[2]*1000/ONE; if(fz<0)fz=-fz;
        printf("[f450-cage] completed. max vertical tracking err %ld mm; ended %ld mm off hover\r\n",
               max_ztrack, fz);
        ok = ok && (max_ztrack < 200) && (fz < 50);   /* tracked the bounce, returned to hover */
    }
    printf(ok ? "[f450-cage] PASS\r\n" : "[f450-cage] FAIL\r\n");
    return ok ? 0 : 1;
}
