/*
 * Cache-scheduled (LPV) TinyMPC for the E1x fabric — Q20 fixed point.
 *
 * Vanilla TinyMPC precomputes one infinite-horizon cache, valid near a single
 * operating point. A fixed-wing aircraft flies across a whole envelope (airspeed,
 * climb, bank), so this example carries a *grid* of caches — one per trimmed
 * flight condition (here: bank × climb coordinated turns of an Ultra Stick 25e) —
 * and at runtime selects the nearest one for the live condition, then runs the
 * same ADMM solve on it. Scheduling is just an array index, so the online cost is
 * identical to single-LTI TinyMPC; the only price is the read-only cache LUT.
 *
 * 8 states [u v w, p q r, phi theta], 4 controls [throttle, elev, ail, rud],
 * horizon 15, 15 ADMM iters, all int32 Q20.
 */
#include <eff.h>
#include "tinympc_lpv.h"

#define NX UT_NX
#define NU UT_NU
#define SH UT_SCALE

int lpv_select_node(int32_t bank_q20, int32_t climb_q20)
{
    int best = 0; int64_t bd = (int64_t)1 << 62;
    for (int n = 0; n < UT_NODES; n++) {
        int64_t db = (int64_t)bank_q20  - UT_node_bank_climb[n][0];
        int64_t dc = (int64_t)climb_q20 - UT_node_bank_climb[n][1];
        int64_t dist = db*db + 9*dc*dc;   /* weight climb vs bank; units are Q20 */
        if (dist < bd) { bd = dist; best = n; }
    }
    return best;
}

__efficient__
void lpv_solve_fabric(const int32_t *e0,
    const int32_t *A, const int32_t *B, const int32_t *Kinf,
    const int32_t *Quu, const int32_t *AmBKt,
    const int32_t *umin, const int32_t *umax, int32_t rho, int N, int iters,
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
            __effcc_parallel(1) for (int i = 0; i < NX; i++) {
                int64_t s=0; for (int m=0;m<NX;m++) s+=(int64_t)AmBKt[i*NX+m]*p[(k+1)*NX+m];
                for (int m=0;m<NU;m++) s-=(int64_t)Kinf[m*NX+i]*r[m];
                p[k*NX+i]=(int32_t)(s>>SH);
            }
        }
        for (int i = 0; i < NX; i++) x[i] = e0[i];
        for (int k = 0; k < N; k++) {
            for (int i = 0; i < NU; i++) { int64_t s=-(int64_t)d[k*NU+i]<<SH; for (int m=0;m<NX;m++) s-=(int64_t)Kinf[i*NX+m]*x[k*NX+m]; u[k*NU+i]=(int32_t)(s>>SH); }
            __effcc_parallel(1) for (int i = 0; i < NX; i++) {
                int64_t s=0; for (int m=0;m<NX;m++) s+=(int64_t)A[i*NX+m]*x[k*NX+m];
                for (int m=0;m<NU;m++) s+=(int64_t)B[i*NU+m]*u[k*NU+m];
                x[(k+1)*NX+i]=(int32_t)(s>>SH);
            }
        }
        for (int k = 0; k < N; k++) for (int i = 0; i < NU; i++) {
            int32_t v=u[k*NU+i]+y[k*NU+i];
            int32_t zn=v<umin[i]?umin[i]:(v>umax[i]?umax[i]:v);
            z[k*NU+i]=zn; y[k*NU+i]+=u[k*NU+i]-zn;
        }
    }
}

/* Scalar reference — identical Q20 arithmetic, no fabric parallelism. */
void lpv_solve_scalar(const int32_t *e0,
    const int32_t *A, const int32_t *B, const int32_t *Kinf,
    const int32_t *Quu, const int32_t *AmBKt,
    const int32_t *umin, const int32_t *umax, int32_t rho, int N, int iters,
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
            int32_t zn=v<umin[i]?umin[i]:(v>umax[i]?umax[i]:v);
            z[k*NU+i]=zn; y[k*NU+i]+=u[k*NU+i]-zn;
        }
    }
}
