#include <stdint.h>
#include <stdlib.h>

#include "dmv_fp.h"

#ifdef EFF_BLD_HAND_OPTIMIZED
__efficient__ void dmv_fp(float_type *a, float_type *b, float_type *z,
                          uint32_t m, uint32_t n, uint32_t mstride)
{
    for (uint32_t i = 0; i < m; i += 8)
    {
        float_type w0 = 0, w1 = 0, w2 = 0, w3 = 0, w4 = 0, w5 = 0, w6 = 0,
                   w7 = 0;
        float_type *a0 = a + mstride * i;
        float_type *a1 = a0 + mstride;
        float_type *a2 = a1 + mstride;
        float_type *a3 = a2 + mstride;
        float_type *a4 = a3 + mstride;
        float_type *a5 = a4 + mstride;
        float_type *a6 = a5 + mstride;
        float_type *a7 = a6 + mstride;
        for (uint32_t j = 0; j < n; j++)
        {
            // Reuse element of b to reduce number of loads in inner loop.
            // Allows all inner loop loads to be placed in NUPEA domain 0.
            float_type v = b[j];
            w0 += a0[j] * v;
            w1 += a1[j] * v;
            w2 += a2[j] * v;
            w3 += a3[j] * v;
            w4 += a4[j] * v;
            w5 += a5[j] * v;
            w6 += a6[j] * v;
            w7 += a7[j] * v;
        }
        z[i] = w0;
        z[i + 1] = w1;
        z[i + 2] = w2;
        z[i + 3] = w3;
        z[i + 4] = w4;
        z[i + 5] = w5;
        z[i + 6] = w6;
        z[i + 7] = w7;
    }
}
#else
__efficient__ void dmv_fp(float_type *a, float_type *b, float_type *z,
                          uint32_t m, uint32_t n, uint32_t mstride)
{
    for (uint32_t i = 0; i < m; ++i)
    {
        float_type w = 0;
        for (uint32_t j = 0; j < n; j++)
        {
            w += a[i * mstride + j] * b[j];
        }
        z[i] = w;
    }
}

#endif
