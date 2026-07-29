#ifdef __EFFCC__
#include <eff.h>
#endif

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "cholesky.h"

__attribute__((always_inline)) int32_t isqrt(int32_t x)
{
    int32_t q = 1, r = 0;
    while (q <= x)
    {
        q <<= 2;
    }
    while (q > 1)
    {
        int32_t t;
        q >>= 2;
        t = x - r - q;
        r >>= 1;
        if (t >= 0)
        {
            x = t;
            r += q;
        }
    }
    return r;
}

#ifdef EFF_BLD_HAND_OPTIMIZED
__efficient__ void cholesky_helper(int *orig, int n, int *chol, int ofs,
                                   int ip1)
{
    __effcc_parallel(5) for (int j = ip1; j < n + ofs; j++)
    {
        int acc = 0;
        for (int k = ofs; k < ip1 - 1; k++)
        {
            acc -= FIXED_POINT_MUL_RESCALE(chol[k * n + ip1 - 1] *
                                           chol[k * n + j]);
        }
        chol[(ip1 - 1) * n + j] = orig[(ip1 - 1) * n + j] + acc;
    }
}

void cholesky(int *orig, int n, int *chol, int ofs)
{
    for (int i = ofs; i < n + ofs; i++)
    {
        chol[i * n + i] = orig[i * n + i];
        for (int k = ofs; k < i; k++)
        {
            chol[i * n + i] -=
                FIXED_POINT_MUL_RESCALE(chol[k * n + i] * chol[k * n + i]);
        }
        chol[i * n + i] = FIXED_POINT_DIVIDEND(isqrt(chol[i * n + i])) >>
                          (FIXED_POINT_BITS / 2);

        cholesky_helper(orig, n, chol, ofs, i + 1);

        for (int j = i + 1; j < n + ofs; j++)
        {
            chol[i * n + j] =
                FIXED_POINT_DIVIDEND(chol[i * n + j]) / chol[i * n + i];
        }
    }
}
#else
__efficient__ void cholesky(int *orig, int n, int *chol, int ofs)
{
    for (int i = ofs; i < n + ofs; i++)
    {
        chol[i * n + i] = orig[i * n + i];
        for (int k = ofs; k < i; k++)
        {
            chol[i * n + i] -=
                FIXED_POINT_MUL_RESCALE(chol[k * n + i] * chol[k * n + i]);
        }
        chol[i * n + i] = FIXED_POINT_DIVIDEND(isqrt(chol[i * n + i])) >>
                          (FIXED_POINT_BITS / 2);

        for (int j = i + 1; j < n + ofs; j++)
        {
            chol[i * n + j] = orig[i * n + j];
            for (int k = ofs; k < i; k++)
            {
                chol[i * n + j] -=
                    FIXED_POINT_MUL_RESCALE(chol[k * n + i] * chol[k * n + j]);
            }
            chol[i * n + j] =
                FIXED_POINT_DIVIDEND(chol[i * n + j]) / chol[i * n + i];
        }
    }
}
#endif
