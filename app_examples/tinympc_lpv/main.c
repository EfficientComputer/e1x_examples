/*
 * Cache-scheduled (LPV) TinyMPC — fixed-wing MPC across a flight envelope.
 *
 * Picks the cache nearest a commanded flight condition (here a 30° right
 * coordinated turn), solves one ADMM step on the fabric and on a bit-identical
 * scalar reference from a small attitude error, and verifies they match. The box
 * is the actuator range minus the selected node's trim control (so the solve
 * optimizes the control *delta* about that trim).
 */
#include <eff.h>
#include <stdio.h>
#include <stdint.h>
#include "tinympc_lpv.h"

#define ONE (1 << UT_SCALE)
#define DEG ((int32_t)(ONE * 0.0174532925))    /* 1 degree in Q20 radians */

/* actuator absolute box (Q20): throttle [0,1], surfaces ±25° */
#define SURF ((int32_t)(ONE * 0.43633))
static const int32_t ABS_MIN[UT_NU] = {0, -SURF, -SURF, -SURF};
static const int32_t ABS_MAX[UT_NU] = {ONE, SURF, SURF, SURF};

static int32_t z0[UT_HOR*UT_NU], y0[UT_HOR*UT_NU], p0[(UT_HOR+1)*UT_NX], d0[UT_HOR*UT_NU], x0[(UT_HOR+1)*UT_NX], u0[UT_HOR*UT_NU];
static int32_t z1[UT_HOR*UT_NU], y1[UT_HOR*UT_NU], p1[(UT_HOR+1)*UT_NX], d1[UT_HOR*UT_NU], x1[(UT_HOR+1)*UT_NX], u1[UT_HOR*UT_NU];

int main(void)
{
    /* commanded flight condition: 30° right bank, level (climb 0) */
    int32_t bank = 30*DEG, climb = 0;
    int node = lpv_select_node(bank, climb);

    /* per-control delta box about this node's trim */
    int32_t umin[UT_NU], umax[UT_NU];
    for (int i = 0; i < UT_NU; i++) { umin[i] = ABS_MIN[i]-UT_utrim[node][i]; umax[i] = ABS_MAX[i]-UT_utrim[node][i]; }

    /* small attitude error: a few deg of roll + pitch off the turn trim */
    int32_t e0[UT_NX] = {0,0,0, 0,0,0, 5*DEG, 2*DEG};

    lpv_solve_fabric(e0, UT_A[node], UT_B[node], UT_Kinf[node], UT_Quu[node], UT_AmBKt[node],
        umin, umax, UT_RHO, UT_HOR, UT_ITERS, z0,y0,p0,d0,x0,u0);
    lpv_solve_scalar(e0, UT_A[node], UT_B[node], UT_Kinf[node], UT_Quu[node], UT_AmBKt[node],
        umin, umax, UT_RHO, UT_HOR, UT_ITERS, z1,y1,p1,d1,x1,u1);

    printf("[tinympc_lpv] Ultra Stick fixed-wing: %d states, %d controls, N=%d, %d turn caches\r\n",
           UT_NX, UT_NU, UT_HOR, UT_NODES);
    printf("[tinympc_lpv] scheduled node %d for bank 30 deg, climb 0\r\n", node);
    printf("[tinympc_lpv] control delta u[0] (x1e6): ");
    for (int i = 0; i < UT_NU; i++) printf("%ld ", (long)z0[i] * 1000000 / ONE);
    printf("\r\n");

    int ok = 1;
    for (int i = 0; i < UT_HOR*UT_NU; i++) if (z0[i] != z1[i]) { ok = 0;
        printf("[tinympc_lpv] MISMATCH at z[%d]: fabric %ld != scalar %ld\r\n", i, (long)z0[i], (long)z1[i]); break; }
    printf(ok ? "[tinympc_lpv] PASS\r\n" : "[tinympc_lpv] FAIL\r\n");
    return ok ? 0 : 1;
}
