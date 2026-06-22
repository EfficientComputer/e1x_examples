#ifndef TINYMPC_LPV_H
#define TINYMPC_LPV_H
#include <stdint.h>
#include "ustick_turn_q20.h"   /* UT_NX, UT_NU, UT_NODES, UT_* LUT */

/* Pick the cache node nearest a flight condition (bank, climb angle), both Q20
 * radians. Trivial vs the solve — this is the "scheduling" of LPV-MPC. */
int  lpv_select_node(int32_t bank_q20, int32_t climb_q20);

/* One scheduled TinyMPC ADMM solve (Q20) on the chosen node's cache; per-control
 * delta box umin/umax. First control delta returned in z[0..UT_NU-1]. */
void lpv_solve_fabric(const int32_t *e0,
    const int32_t *A, const int32_t *B, const int32_t *Kinf,
    const int32_t *Quu, const int32_t *AmBKt,
    const int32_t *umin, const int32_t *umax, int32_t rho, int N, int iters,
    int32_t *z, int32_t *y, int32_t *p, int32_t *d, int32_t *x, int32_t *u);

void lpv_solve_scalar(const int32_t *e0,
    const int32_t *A, const int32_t *B, const int32_t *Kinf,
    const int32_t *Quu, const int32_t *AmBKt,
    const int32_t *umin, const int32_t *umax, int32_t rho, int N, int iters,
    int32_t *z, int32_t *y, int32_t *p, int32_t *d, int32_t *x, int32_t *u);

#endif
