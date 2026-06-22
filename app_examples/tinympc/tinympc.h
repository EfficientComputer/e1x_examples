#ifndef TINYMPC_H
#define TINYMPC_H
#include <stdint.h>

/* Crazyflie-class quadrotor MPC: 12 states, 4 controls, horizon N=10. */
#define NX 12
#define NU 4
#define NHORIZON 10
#define RHO 5            /* integer ADMM penalty (kept in Q20 by the recursion) */

/* One TinyMPC ADMM solve (all iterations) in Q20 fixed point. e0 is the initial
 * error state (x - xref) in Q20; the first control is returned in z[0..NU-1].
 * Caller provides the cache (A,B,Kinf,Quu,AmBKt) and scratch (z,y,p,d,x,u). */
void tinympc_solve_fabric(const int32_t *e0,
    const int32_t *A, const int32_t *B, const int32_t *Kinf,
    const int32_t *Quu, const int32_t *AmBKt, int32_t umin, int32_t umax,
    int32_t rho, int N, int iters,
    int32_t *z, int32_t *y, int32_t *p, int32_t *d, int32_t *x, int32_t *u);

/* Bit-identical scalar reference (same Q20 math, no fabric parallelism). */
void tinympc_solve_scalar(const int32_t *e0,
    const int32_t *A, const int32_t *B, const int32_t *Kinf,
    const int32_t *Quu, const int32_t *AmBKt, int32_t umin, int32_t umax,
    int32_t rho, int N, int iters,
    int32_t *z, int32_t *y, int32_t *p, int32_t *d, int32_t *x, int32_t *u);

#endif
