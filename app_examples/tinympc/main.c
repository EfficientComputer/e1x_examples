/*
 * TinyMPC quadrotor MPC — E1x demo + correctness check.
 *
 * Runs one TinyMPC ADMM solve on the official Crazyflie-class quadrotor cache
 * (12 states, 4 controls, horizon 10) from a hover error state, on the fabric
 * kernel and on a bit-identical scalar reference, and verifies they match.
 * The first control z[0..3] is the command the controller would apply this step.
 */
#include <eff.h>
#include <stdio.h>
#include <stdint.h>
#include "tinympc.h"
#include "quad_cache_q20.h"   /* QC_*_Q20 cache (Q20), QCQ=20, QC_UMIN/UMAX_Q20 */

#define ONE (1 << QCQ)

/* scratch (one set per solver so we can compare) */
static int32_t z0[NHORIZON*NU], y0[NHORIZON*NU], p0[(NHORIZON+1)*NX], d0[NHORIZON*NU], x0[(NHORIZON+1)*NX], u0[NHORIZON*NU];
static int32_t z1[NHORIZON*NU], y1[NHORIZON*NU], p1[(NHORIZON+1)*NX], d1[NHORIZON*NU], x1[(NHORIZON+1)*NX], u1[NHORIZON*NU];

int main(void)
{
    /* hover error state e0 = x - xref (Q20): 1 m in y, plus small attitude/vel */
    int32_t e0[NX] = {0, ONE, 0, ONE/5, 0, 0, ONE/10, 0, 0, 0, 0, 0};

    /* fabric solve */
    tinympc_solve_fabric(e0, QC_A_Q20, QC_B_Q20, QC_Kinf_Q20, QC_Quu_Q20, QC_AmBKt_Q20,
        QC_UMIN_Q20, QC_UMAX_Q20, RHO, NHORIZON, NHORIZON, z0, y0, p0, d0, x0, u0);
    /* scalar reference */
    tinympc_solve_scalar(e0, QC_A_Q20, QC_B_Q20, QC_Kinf_Q20, QC_Quu_Q20, QC_AmBKt_Q20,
        QC_UMIN_Q20, QC_UMAX_Q20, RHO, NHORIZON, NHORIZON, z1, y1, p1, d1, x1, u1);

    printf("[tinympc] Crazyflie quad MPC: %d states, %d controls, N=%d, %d ADMM iters (Q20)\r\n",
           NX, NU, NHORIZON, NHORIZON);
    printf("[tinympc] control u[0] (x1e6): ");
    for (int i = 0; i < NU; i++) {
        long micro = (long)z0[i] * 1000000 / ONE;   /* Q20 -> millionths, integer */
        printf("%ld ", micro);
    }
    printf("\r\n");

    /* correctness: fabric must match the scalar reference exactly */
    int ok = 1;
    for (int i = 0; i < NHORIZON*NU; i++) if (z0[i] != z1[i]) { ok = 0;
        printf("[tinympc] MISMATCH at z[%d]: fabric %ld != scalar %ld\r\n", i, (long)z0[i], (long)z1[i]); break; }
    printf(ok ? "[tinympc] PASS\r\n" : "[tinympc] FAIL\r\n");
    return ok ? 0 : 1;
}
