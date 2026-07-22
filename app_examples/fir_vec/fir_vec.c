#ifdef __EFFCC__
#include <eff.h>
#endif

#include <stdint.h>

#include "fir_vec.h"

// vdot_s16 vectorized FIR: pack 16 taps into 8 _v2s16_t pairs.
// Per output sample: 8 vdot_s16 calls vs 16 scalar MACs.
// Assumes W=16, values fit in int16.
__efficient__ void fir(int *x, int *w, int *o)
{
    // Pre-pack filter weights into 8 pairs — loaded once per call
    _v2s16_t wv0 = {(int16_t)w[0], (int16_t)w[1]};
    _v2s16_t wv1 = {(int16_t)w[2], (int16_t)w[3]};
    _v2s16_t wv2 = {(int16_t)w[4], (int16_t)w[5]};
    _v2s16_t wv3 = {(int16_t)w[6], (int16_t)w[7]};
    _v2s16_t wv4 = {(int16_t)w[8], (int16_t)w[9]};
    _v2s16_t wv5 = {(int16_t)w[10], (int16_t)w[11]};
    _v2s16_t wv6 = {(int16_t)w[12], (int16_t)w[13]};
    _v2s16_t wv7 = {(int16_t)w[14], (int16_t)w[15]};

    // Prime the shift register with the first 15 input samples
    int r0 = x[0], r1 = x[1], r2 = x[2], r3 = x[3];
    int r4 = x[4], r5 = x[5], r6 = x[6], r7 = x[7];
    int r8 = x[8], r9 = x[9], r10 = x[10], r11 = x[11];
    int r12 = x[12], r13 = x[13], r14 = x[14];

    for (int i = 0; i < N; ++i)
    {
        int r15 = x[i + 15];

        _v2s16_t rv0 = {(int16_t)r0, (int16_t)r1};
        _v2s16_t rv1 = {(int16_t)r2, (int16_t)r3};
        _v2s16_t rv2 = {(int16_t)r4, (int16_t)r5};
        _v2s16_t rv3 = {(int16_t)r6, (int16_t)r7};
        _v2s16_t rv4 = {(int16_t)r8, (int16_t)r9};
        _v2s16_t rv5 = {(int16_t)r10, (int16_t)r11};
        _v2s16_t rv6 = {(int16_t)r12, (int16_t)r13};
        _v2s16_t rv7 = {(int16_t)r14, (int16_t)r15};

        o[i] = __builtin_effcc_vdot_s16(wv0, rv0) +
               __builtin_effcc_vdot_s16(wv1, rv1) +
               __builtin_effcc_vdot_s16(wv2, rv2) +
               __builtin_effcc_vdot_s16(wv3, rv3) +
               __builtin_effcc_vdot_s16(wv4, rv4) +
               __builtin_effcc_vdot_s16(wv5, rv5) +
               __builtin_effcc_vdot_s16(wv6, rv6) +
               __builtin_effcc_vdot_s16(wv7, rv7);

        r0 = r1;
        r1 = r2;
        r2 = r3;
        r3 = r4;
        r4 = r5;
        r5 = r6;
        r6 = r7;
        r7 = r8;
        r8 = r9;
        r9 = r10;
        r10 = r11;
        r11 = r12;
        r12 = r13;
        r13 = r14;
        r14 = r15;
    }
}
