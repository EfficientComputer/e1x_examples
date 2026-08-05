#include "mekf_q20.h"
#include <string.h>
#include <math.h>

/* ── Multiply helpers ─────────────────────────────────────────────────── */

/* Q15 × Q15 → Q15 */
static inline int32_t mul15(int32_t a, int32_t b) { return (a * b) >> 15; }

/* Q15 × Q20 → Q20  (dt_w × P element; |a|≤164 so no overflow) */
static inline int32_t mul15x20(int32_t a_q15, int32_t b_q20)
{
    return (a_q15 * b_q20) >> 15;
}

/* Q8 × Q20 → Q20  (h element × P element; |a|≤256, |b|≤94372 → 24M < 2^31) */
static inline int32_t mul8x20(int32_t a_q8, int32_t b_q20)
{
    return (a_q8 * b_q20) >> 8;
}

/* Q15 × Q13 → Q20  (PHT_q15 × Sinv_q13, then >>8; used in K computation)
 * |a|≤2949, |b|≤87000 → 256M < 2^31 ✓ */
static inline int32_t mul15x13_to_q20(int32_t a_q15, int32_t b_q13)
{
    return (a_q15 * b_q13) >> 8;   /* Q15+Q13=Q28, >>8 → Q20 */
}

/* Q8 × Q20 → Q20  (K_q8 × PHT_q20 for P update; |a|≤256, |b|≤94372) */
static inline int32_t mul8x20_p(int32_t a_q8, int32_t b_q20)
{
    return (a_q8 * b_q20) >> 8;
}

/* ── Quaternion helpers (Q15) ─────────────────────────────────────────── */

static void q_norm(int32_t q[4])
{
    int32_t sq = 0;
    for (int i = 0; i < 4; i++) sq += (q[i] * q[i]) >> 15;
    int32_t half_err = (sq - Q15_ONE) >> 1;
    for (int i = 0; i < 4; i++) q[i] -= mul15(q[i], half_err);
}

static void q_mul(int32_t out[4], const int32_t p[4], const int32_t r[4])
{
    out[0] = mul15(p[0],r[0]) - mul15(p[1],r[1]) - mul15(p[2],r[2]) - mul15(p[3],r[3]);
    out[1] = mul15(p[0],r[1]) + mul15(p[1],r[0]) + mul15(p[2],r[3]) - mul15(p[3],r[2]);
    out[2] = mul15(p[0],r[2]) - mul15(p[1],r[3]) + mul15(p[2],r[0]) + mul15(p[3],r[1]);
    out[3] = mul15(p[0],r[3]) + mul15(p[1],r[2]) - mul15(p[2],r[1]) + mul15(p[3],r[0]);
}

/* ── 3×3 inverse in Q15 (Cramer's rule) ───────────────────────────────── */
static void inv3_q15(int32_t Sinv[9], const int32_t S[9])
{
    int32_t c00 = mul15(S[4],S[8]) - mul15(S[5],S[7]);
    int32_t c01 = mul15(S[5],S[6]) - mul15(S[3],S[8]);
    int32_t c02 = mul15(S[3],S[7]) - mul15(S[4],S[6]);
    int32_t c10 = mul15(S[2],S[7]) - mul15(S[1],S[8]);
    int32_t c11 = mul15(S[0],S[8]) - mul15(S[2],S[6]);
    int32_t c12 = mul15(S[1],S[6]) - mul15(S[0],S[7]);
    int32_t c20 = mul15(S[1],S[5]) - mul15(S[2],S[4]);
    int32_t c21 = mul15(S[2],S[3]) - mul15(S[0],S[5]);
    int32_t c22 = mul15(S[0],S[4]) - mul15(S[1],S[3]);
    int32_t det = mul15(S[0],c00) + mul15(S[1],c01) + mul15(S[2],c02);
    if (det == 0) det = 1;
    Sinv[0]=(c00<<15)/det; Sinv[1]=(c10<<15)/det; Sinv[2]=(c20<<15)/det;
    Sinv[3]=(c01<<15)/det; Sinv[4]=(c11<<15)/det; Sinv[5]=(c21<<15)/det;
    Sinv[6]=(c02<<15)/det; Sinv[7]=(c12<<15)/det; Sinv[8]=(c22<<15)/det;
}

/* ── Init ──────────────────────────────────────────────────────────────── */
void mekf_q20_init(MEKF_Q20 *m,
                   const float q0[4], const float bg0[3],
                   float sigma_q0, float sigma_bg0)
{
    for (int i = 0; i < 4; i++) m->q[i]   = (int32_t)(q0[i]   * Q15_ONE);
    for (int i = 0; i < 3; i++) m->b_g[i] = (int32_t)(bg0[i]  * Q15_ONE);
    for (int i = 0; i < 36; i++) m->P[i]  = 0;
    int32_t sq0  = (int32_t)(sigma_q0  * sigma_q0  * Q20_ONE);
    int32_t sbg0 = (int32_t)(sigma_bg0 * sigma_bg0 * Q20_ONE);
    m->P[0*6+0]=sq0;  m->P[1*6+1]=sq0;  m->P[2*6+2]=sq0;
    m->P[3*6+3]=sbg0; m->P[4*6+4]=sbg0; m->P[5*6+5]=sbg0;
}

/* ── Predict — P in Q20, dt/w in Q15 ──────────────────────────────────── */
void mekf_q20_predict(MEKF_Q20 *m,
                      const int32_t w_gyro_q15[3],
                      int32_t dt_q15,
                      int32_t qg_q20,
                      int32_t qbg_q20)
{
    int32_t w[3];
    w[0] = w_gyro_q15[0] - m->b_g[0];
    w[1] = w_gyro_q15[1] - m->b_g[1];
    w[2] = w_gyro_q15[2] - m->b_g[2];

    /* Quaternion kinematics: dq = [1, w*dt/2], Q15.
     * Round-to-nearest via (x + Q15_ONE) >> 16 avoids the systematic
     * negative-drift from arithmetic right-shift rounding toward -inf. */
    int32_t dq[4];
    dq[0] = Q15_ONE;
    dq[1] = (w[0] * dt_q15 + Q15_ONE) >> 16;
    dq[2] = (w[1] * dt_q15 + Q15_ONE) >> 16;
    dq[3] = (w[2] * dt_q15 + Q15_ONE) >> 16;

    int32_t q_new[4];
    q_mul(q_new, m->q, dq);
    q_norm(q_new);
    m->q[0]=q_new[0]; m->q[1]=q_new[1]; m->q[2]=q_new[2]; m->q[3]=q_new[3];

    /* dt*w terms in Q15 for Phi computation */
    int32_t dt_w0 = mul15(dt_q15, w[0]);
    int32_t dt_w1 = mul15(dt_q15, w[1]);
    int32_t dt_w2 = mul15(dt_q15, w[2]);

    /* ΦP in Q20: top rows = (I−[w×]dt)·P[0:3,:]−dt·P[3:6,:]
     * mul15x20: Q15×Q20→Q20  (|dt_w|≤164, |P|≤94372 → 15.5M < 2^31) */
    int32_t PhiP[36];
    for (int j = 0; j < 6; j++) {
        int32_t p0=m->P[0*6+j], p1=m->P[1*6+j], p2=m->P[2*6+j];
        int32_t p3=m->P[3*6+j], p4=m->P[4*6+j], p5=m->P[5*6+j];
        PhiP[3*6+j] = p3;
        PhiP[4*6+j] = p4;
        PhiP[5*6+j] = p5;
        PhiP[0*6+j] = p0 + mul15x20(dt_w2,p1) - mul15x20(dt_w1,p2) - mul15x20(dt_q15,p3);
        PhiP[1*6+j] = p1 - mul15x20(dt_w2,p0) + mul15x20(dt_w0,p2) - mul15x20(dt_q15,p4);
        PhiP[2*6+j] = p2 + mul15x20(dt_w1,p0) - mul15x20(dt_w0,p1) - mul15x20(dt_q15,p5);
    }

    /* ΦPΦᵀ in Q20 */
    int32_t P_new[36];
    for (int i = 0; i < 6; i++) {
        int32_t m0=PhiP[i*6+0], m1=PhiP[i*6+1], m2=PhiP[i*6+2];
        int32_t m3=PhiP[i*6+3], m4=PhiP[i*6+4], m5=PhiP[i*6+5];
        P_new[i*6+3] = m3;
        P_new[i*6+4] = m4;
        P_new[i*6+5] = m5;
        P_new[i*6+0] = m0 + mul15x20(dt_w2,m1) - mul15x20(dt_w1,m2) - mul15x20(dt_q15,m3);
        P_new[i*6+1] = m1 - mul15x20(dt_w2,m0) + mul15x20(dt_w0,m2) - mul15x20(dt_q15,m4);
        P_new[i*6+2] = m2 + mul15x20(dt_w1,m0) - mul15x20(dt_w0,m1) - mul15x20(dt_q15,m5);
    }

    memcpy(m->P, P_new, 36*sizeof(int32_t));
    m->P[0*6+0]+=qg_q20;  m->P[1*6+1]+=qg_q20;  m->P[2*6+2]+=qg_q20;
    m->P[3*6+3]+=qbg_q20; m->P[4*6+4]+=qbg_q20; m->P[5*6+5]+=qbg_q20;
}

/* ── Accel update — mixed Q8/Q15/Q20 ──────────────────────────────────── */
void mekf_q20_update_accel(MEKF_Q20 *m,
                            const int32_t a_norm_q8[3],
                            int32_t Ra_q20)
{
    /* h = R^T[:,2] in Q8: third row of rotation matrix × 256.
     * R[2,0]=2(qxqz−qwqy), R[2,1]=2(qyqz+qwqx), R[2,2]=1−2(qx²+qy²)
     * Values in Q15, scale down to Q8 by >>7. */
    int32_t qw=m->q[0], qx=m->q[1], qy=m->q[2], qz=m->q[3];
    int32_t h[3];
    h[0] = (2*(qx*qz - qw*qy)) >> 22;   /* Q30>>22 = Q8 */
    h[1] = (2*(qy*qz + qw*qx)) >> 22;
    h[2] = ((int32_t)Q8_ONE) - ((2*(qx*qx + qy*qy)) >> 22);

    /* Innovation z = a_norm_q8 - h_q8 */
    int32_t z[3] = { a_norm_q8[0]-h[0], a_norm_q8[1]-h[1], a_norm_q8[2]-h[2] };

    /* PHT[6×3] in Q20: P_q20 × H_q8^T, products Q28 >>8 = Q20.
     * H is skew(h): col0=[0,h2,-h1], col1=[-h2,0,h0], col2=[h1,-h0,0]
     * |h_q8|≤256, |P_q20|≤94372 → 24M < 2^31 ✓ */
    int32_t PHT[18];
    for (int i = 0; i < 6; i++) {
        int32_t p0=m->P[i*6+0], p1=m->P[i*6+1], p2=m->P[i*6+2];
        PHT[i*3+0] = mul8x20(-h[2],p1) + mul8x20(h[1],p2);
        PHT[i*3+1] = mul8x20( h[2],p0) + mul8x20(-h[0],p2);
        PHT[i*3+2] = mul8x20(-h[1],p0) + mul8x20(h[0],p1);
    }

    /* S_q20[3×3] = H_q8 × PHT_q20[0:3,:] + Ra_q20·I
     * Same overflow budget as PHT: 24M < 2^31 ✓ */
    int32_t S_q20[9];
    for (int j = 0; j < 3; j++) {
        S_q20[0*3+j] = mul8x20(-h[2],PHT[1*3+j]) + mul8x20(h[1],PHT[2*3+j]);
        S_q20[1*3+j] = mul8x20( h[2],PHT[0*3+j]) + mul8x20(-h[0],PHT[2*3+j]);
        S_q20[2*3+j] = mul8x20(-h[1],PHT[0*3+j]) + mul8x20(h[0],PHT[1*3+j]);
    }
    S_q20[0]+=Ra_q20; S_q20[4]+=Ra_q20; S_q20[8]+=Ra_q20;

    /* Regularize: add h_q8 × h_q8^T to inflate the unobservable
     * rotation-about-gravity eigenvalue (same as float MEKF).
     * h_q8[i]*h_q8[j] is Q16; <<4 converts to Q20. */
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            S_q20[i*3+j] += (h[i] * h[j]) << 4;

    /* Diagonal-only Sinv: Sinv_q0[j] = Q20_ONE / S_q20[j,j].
     * Avoids inv3 entirely. Valid since S is nearly diagonal for MEKF.
     * Sinv_real[j] = 1/S_real[j] ≈ 10.5; in Q0 (integer) ≈ 10-160.
     * Division is exact int32, no overflow: Q20_ONE=1048576 < 2^31. */
    int32_t Sinv_q0[3];
    {
        int32_t s0=S_q20[0], s4=S_q20[4], s8=S_q20[8];
        Sinv_q0[0] = (s0 > 0) ? Q20_ONE / s0 : 1;
        Sinv_q0[1] = (s4 > 0) ? Q20_ONE / s4 : 1;
        Sinv_q0[2] = (s8 > 0) ? Q20_ONE / s8 : 1;
    }

    /* K_q15[6×3] = (PHT_q20 × Sinv_q0) >> 5.
     * PHT_q20 × Sinv_q0 = K_real × Q20_ONE ≤ Q20_ONE (since K ≤ 1).
     * >>5 gives K in Q15; max = Q15_ONE = 32768 ✓
     * No overflow: PHT×Sinv ≤ Q20_ONE = 1048576 (bounded since K ≤ 1) */
    int32_t K_q15[18];
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 3; j++)
            K_q15[i*3+j] = (PHT[i*3+j] * Sinv_q0[j]) >> 5;

    /* K projection: remove rotation-about-gravity (unobservable yaw) from
     * the attitude rows. Same as float MEKF's explicit projection.
     * dot_q15 = h^T × K_att[:,j] (Q8×Q15=Q23, >>8 → Q15).
     * K_att[i,j] -= h[i] × dot_q15 >> 8 (Q8×Q15=Q23, >>8 → Q15). */
    for (int j = 0; j < 3; j++) {
        int32_t dot = ((int32_t)h[0]*K_q15[0*3+j] + (int32_t)h[1]*K_q15[1*3+j]
                     + (int32_t)h[2]*K_q15[2*3+j]) >> 8;
        for (int i = 0; i < 3; i++)
            K_q15[i*3+j] -= (h[i] * dot) >> 8;
    }

    /* dx_q20 = K_q15 × z_q8 >> 3.  Q15×Q8=Q23; >>3 → Q20.
     * Sum of 3: 3×32768×256=25M < 2^31 ✓ */
    int32_t dx_q20[6] = {0};
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 3; j++)
            dx_q20[i] += (K_q15[i*3+j] * z[j]) >> 3;

    /* Attitude correction: dq = [1, dx/2] in Q15.
     * dx_q20 is Q20; >>5 = Q15; >>1 for half-angle = >>6 total. */
    int32_t dq_c[4] = { Q15_ONE, dx_q20[0]>>5>>1, dx_q20[1]>>5>>1, dx_q20[2]>>5>>1 };
    int32_t q_new[4];
    q_mul(q_new, m->q, dq_c);
    q_norm(q_new);
    m->q[0]=q_new[0]; m->q[1]=q_new[1]; m->q[2]=q_new[2]; m->q[3]=q_new[3];

    /* Bias correction: dx_q20>>5 = dx_q15 added to b_g_q15 */
    m->b_g[0] += dx_q20[3] >> 5;
    m->b_g[1] += dx_q20[4] >> 5;
    m->b_g[2] += dx_q20[5] >> 5;

    /* P update: P_q20 -= K × PHT^T for att rows (i<3).
     * K_q20≤Q20_ONE, PHT_q20≤Q20_ONE → product Q40 may overflow.
     * Use K_q15 = K_q20>>5, PHT_q15 = PHT_q20>>5, result >>10 = Q20.
     * |K_q15|≤32768, |PHT_q15|≤Q20_ONE>>5 → product ≤ 32768×Q20_ONE/32 < 2^31 ✓ */
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 6; j++) {
            int32_t s = 0;
            for (int k = 0; k < 3; k++)
                s += ((K_q15[i*3+k]) * (PHT[j*3+k] >> 5)) >> 10;
            m->P[i*6+j] -= s;
        }
    }
    /* Mirror att-row updates to bias cross-terms */
    for (int i = 0; i < 3; i++)
        for (int j = 3; j < 6; j++)
            m->P[j*6+i] = m->P[i*6+j];
    /* Symmetrize att block */
    for (int i = 0; i < 3; i++)
        for (int j = i+1; j < 3; j++) {
            int32_t avg = (m->P[i*6+j] + m->P[j*6+i]) >> 1;
            m->P[i*6+j] = m->P[j*6+i] = avg;
        }
    /* Clamp diagonal to ≥ 1 to prevent negative variance from rounding */
    for (int i = 0; i < 6; i++)
        if (m->P[i*6+i] < 1) m->P[i*6+i] = 1;
}

/* ── Readback ──────────────────────────────────────────────────────────── */
float mekf_q20_err_deg(const MEKF_Q20 *m, const double q_ref[4])
{
    float dot = 0.0f;
    for (int i = 0; i < 4; i++)
        dot += (m->q[i] / (float)Q15_ONE) * (float)q_ref[i];
    if (dot >  1.0f) dot =  1.0f;
    if (dot < -1.0f) dot = -1.0f;
    return 2.0f * acosf(fabsf(dot)) * 180.0f / 3.14159265f;
}

void mekf_q20_get_bg(const MEKF_Q20 *m, float bg[3])
{
    for (int i = 0; i < 3; i++) bg[i] = m->b_g[i] / (float)Q15_ONE;
}
