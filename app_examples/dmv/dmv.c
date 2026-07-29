#ifdef __EFFCC__
#include <eff.h>
#endif

#include "dmv.h"

#include <stdint.h>
#include <stdlib.h>

__efficient__ void dmv(int32_t *a, int32_t *b, int32_t *z, uint32_t n,
                       uint32_t m)
{
#ifdef EFF_BLD_HAND_OPTIMIZED
    for (uint32_t i = 0; i < n; i += 4)
    {
        int32_t w0 = 0, w1 = 0, w2 = 0, w3 = 0;
        uint32_t row_idx = n * i;
        int32_t *a0 = a + row_idx;
        int32_t *a1 = a0 + n;
        int32_t *a2 = a1 + n;
        int32_t *a3 = a2 + n;
        __effcc_parallel(4) for (uint32_t j = 0; j < m; j++)
        {
            int32_t v = b[j];
            w0 += a0[j] * v;
            w1 += a1[j] * v;
            w2 += a2[j] * v;
            w3 += a3[j] * v;
        }
        z[i] = w0;
        z[i + 1] = w1;
        z[i + 2] = w2;
        z[i + 3] = w3;
    }
#else
    for (uint32_t i = 0; i < n; i++)
    {
        int32_t w = 0;
        for (uint32_t j = 0; j < m; j++)
        {
            w += a[n * i + j] * b[j];
        }
        z[i] = w;
    }
#endif
}
