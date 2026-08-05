/*
 * TinyMPC ADMM solve for the E1x fabric — Q20 fixed point.
 *
 * TinyMPC (https://tinympc.org) is model-predictive control built on a
 * precomputed infinite-horizon cache (Kinf, Pinf, Quu^-1, AmBKt from the DARE)
 * plus an online ADMM loop of constant-matrix matvec sweeps — no online matrix
 * factorization. That structure maps well to the E1x: the whole multi-iteration
 * solve is one __efficient__ fabric dispatch, in int32 Q20 (int32*int32->int64
 * widening accumulate, >>20) with no float and no 64-bit divide.
 *
 * Per ADMM iteration:
 *   backward:  r = -rho*(z - y);  d = Quu^-1 (B^T p + r);  p = AmBKt^T p - Kinf^T r
 *   forward:   u = -Kinf x - d;   x+ = A x + B u
 *   slack+dual: z = clip(u + y, [umin,umax]);  y += u - z
 *
 * This file provides the fabric kernel and a bit-identical scalar reference
 * (same Q20 ops, no parallelism) for the correctness check in main.c.
 */
#include <eff.h>
#include "tinympc.h"

#define SH 20            /* Q20 */

__efficient__
void tinympc_solve_fabric(const int32_t *e0,
    const int32_t *A, const int32_t *B, const int32_t *Kinf,
    const int32_t *Quu, const int32_t *AmBKt, int32_t umin, int32_t umax,
    int32_t rho, int N, int iters,
    int32_t *z, int32_t *y, int32_t *p, int32_t *d, int32_t *x, int32_t *u)
{
    for (int i = 0; i < N*NU; i++) { z[i] = 0; y[i] = 0; }
    for (int it = 0; it < iters; it++) {
        for (int i = 0; i < NX; i++) p[N*NX+i] = 0;
        /* backward */
        for (int k = N-1; k >= 0; k--) {
            int32_t r[NU];
            for (int i = 0; i < NU; i++) r[i] = -rho*(z[k*NU+i]-y[k*NU+i]);
            int32_t Btp[NU];
            for (int i = 0; i < NU; i++) { int64_t s=0; for (int m=0;m<NX;m++) s+=(int64_t)B[m*NU+i]*p[(k+1)*NX+m]; Btp[i]=(int32_t)(s>>SH); }
            for (int i = 0; i < NU; i++) { int64_t s=0; for (int m=0;m<NU;m++) s+=(int64_t)Quu[i*NU+m]*(Btp[m]+r[m]); d[k*NU+i]=(int32_t)(s>>SH); }
            __effcc_parallel(1) for (int i = 0; i < NX; i++) {
                int64_t s=0; for (int m=0;m<NX;m++) s+=(int64_t)AmBKt[i*NX+m]*p[(k+1)*NX+m];
                for (int m=0;m<NU;m++) s-=(int64_t)Kinf[m*NX+i]*r[m];
                p[k*NX+i]=(int32_t)(s>>SH);
            }
        }
        /* forward */
        for (int i = 0; i < NX; i++) x[i] = e0[i];
        for (int k = 0; k < N; k++) {
            for (int i = 0; i < NU; i++) { int64_t s=-(int64_t)d[k*NU+i]<<SH; for (int m=0;m<NX;m++) s-=(int64_t)Kinf[i*NX+m]*x[k*NX+m]; u[k*NU+i]=(int32_t)(s>>SH); }
            __effcc_parallel(1) for (int i = 0; i < NX; i++) {
                int64_t s=0; for (int m=0;m<NX;m++) s+=(int64_t)A[i*NX+m]*x[k*NX+m];
                for (int m=0;m<NU;m++) s+=(int64_t)B[i*NU+m]*u[k*NU+m];
                x[(k+1)*NX+i]=(int32_t)(s>>SH);
            }
        }
        /* slack + dual */
        for (int k = 0; k < N; k++) for (int i = 0; i < NU; i++) {
            int32_t v=u[k*NU+i]+y[k*NU+i];
            int32_t zn=v<umin?umin:(v>umax?umax:v);
            z[k*NU+i]=zn; y[k*NU+i]+=u[k*NU+i]-zn;
        }
    }
}

/* Scalar reference — identical arithmetic, no __efficient__ / __effcc_parallel. */
void tinympc_solve_scalar(const int32_t *e0,
    const int32_t *A, const int32_t *B, const int32_t *Kinf,
    const int32_t *Quu, const int32_t *AmBKt, int32_t umin, int32_t umax,
    int32_t rho, int N, int iters,
    int32_t *z, int32_t *y, int32_t *p, int32_t *d, int32_t *x, int32_t *u)
{
    for (int i = 0; i < N*NU; i++) { z[i] = 0; y[i] = 0; }
    for (int it = 0; it < iters; it++) {
        for (int i = 0; i < NX; i++) p[N*NX+i] = 0;
        for (int k = N-1; k >= 0; k--) {
            int32_t r[NU];
            for (int i = 0; i < NU; i++) r[i] = -rho*(z[k*NU+i]-y[k*NU+i]);
            int32_t Btp[NU];
            for (int i = 0; i < NU; i++) { int64_t s=0; for (int m=0;m<NX;m++) s+=(int64_t)B[m*NU+i]*p[(k+1)*NX+m]; Btp[i]=(int32_t)(s>>SH); }
            for (int i = 0; i < NU; i++) { int64_t s=0; for (int m=0;m<NU;m++) s+=(int64_t)Quu[i*NU+m]*(Btp[m]+r[m]); d[k*NU+i]=(int32_t)(s>>SH); }
            for (int i = 0; i < NX; i++) {
                int64_t s=0; for (int m=0;m<NX;m++) s+=(int64_t)AmBKt[i*NX+m]*p[(k+1)*NX+m];
                for (int m=0;m<NU;m++) s-=(int64_t)Kinf[m*NX+i]*r[m];
                p[k*NX+i]=(int32_t)(s>>SH);
            }
        }
        for (int i = 0; i < NX; i++) x[i] = e0[i];
        for (int k = 0; k < N; k++) {
            for (int i = 0; i < NU; i++) { int64_t s=-(int64_t)d[k*NU+i]<<SH; for (int m=0;m<NX;m++) s-=(int64_t)Kinf[i*NX+m]*x[k*NX+m]; u[k*NU+i]=(int32_t)(s>>SH); }
            for (int i = 0; i < NX; i++) {
                int64_t s=0; for (int m=0;m<NX;m++) s+=(int64_t)A[i*NX+m]*x[k*NX+m];
                for (int m=0;m<NU;m++) s+=(int64_t)B[i*NU+m]*u[k*NU+m];
                x[(k+1)*NX+i]=(int32_t)(s>>SH);
            }
        }
        for (int k = 0; k < N; k++) for (int i = 0; i < NU; i++) {
            int32_t v=u[k*NU+i]+y[k*NU+i];
            int32_t zn=v<umin?umin:(v>umax?umax:v);
            z[k*NU+i]=zn; y[k*NU+i]+=u[k*NU+i]-zn;
        }
    }
}
