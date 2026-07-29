#include "fir.h"

__efficient__ void fir(int *x, int *w, int *o)
{
#ifdef EFF_BLD_HAND_OPTIMIZED
    int w0 = w[0], w1 = w[1], w2 = w[2], w3 = w[3];
    int w4 = w[4], w5 = w[5], w6 = w[6], w7 = w[7];
    int w8 = w[8], w9 = w[9], w10 = w[10], w11 = w[11];
    int w12 = w[12], w13 = w[13], w14 = w[14], w15 = w[15];

    int r0 = x[0], r1 = x[1], r2 = x[2], r3 = x[3];
    int r4 = x[4], r5 = x[5], r6 = x[6], r7 = x[7];
    int r8 = x[8], r9 = x[9], r10 = x[10], r11 = x[11];
    int r12 = x[12], r13 = x[13], r14 = x[14];

    for (int i = 0; i < N; ++i)
    {
        int r15 = x[i + 15];
        int s = w0 * r0 + w1 * r1 + w2 * r2 + w3 * r3 + w4 * r4 + w5 * r5 +
                w6 * r6 + w7 * r7 + w8 * r8 + w9 * r9 + w10 * r10 + w11 * r11 +
                w12 * r12 + w13 * r13 + w14 * r14 + w15 * r15;
        o[i] = s;
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
#else
    for (int i = 0; i < N; ++i) {
        int sum = 0;
        for (int j = 0; j < W; ++j) {
            sum += w[j] * x[i + j];
        }
        o[i] = sum;
    }
#endif
}
