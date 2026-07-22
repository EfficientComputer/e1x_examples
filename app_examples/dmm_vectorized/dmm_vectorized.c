#ifdef __EFFCC__
#include <effintrinsic.h>
#endif

#include <stdint.h>
#include <stdlib.h>

#include "dmm_vectorzied.h"

// 4x4 block transpose: reads 4 sequential elements per B row, writes 4
// sequential elements per B_T row — reduces distinct router paths vs scalar
// element-wise.
__efficient__ void transpose_b(const int16_t *restrict b, int16_t *restrict b_T,
                               uint32_t k, uint32_t n)
{
    for (uint32_t j = 0; j < n; j += 4)
    {
        __effcc_parallel(2) for (uint32_t i = 0; i < k; i += 4)
        {
            int16_t t00 = b[(i + 0) * n + j + 0], t01 = b[(i + 0) * n + j + 1],
                    t02 = b[(i + 0) * n + j + 2], t03 = b[(i + 0) * n + j + 3];
            int16_t t10 = b[(i + 1) * n + j + 0], t11 = b[(i + 1) * n + j + 1],
                    t12 = b[(i + 1) * n + j + 2], t13 = b[(i + 1) * n + j + 3];
            int16_t t20 = b[(i + 2) * n + j + 0], t21 = b[(i + 2) * n + j + 1],
                    t22 = b[(i + 2) * n + j + 2], t23 = b[(i + 2) * n + j + 3];
            int16_t t30 = b[(i + 3) * n + j + 0], t31 = b[(i + 3) * n + j + 1],
                    t32 = b[(i + 3) * n + j + 2], t33 = b[(i + 3) * n + j + 3];
            b_T[(j + 0) * k + i + 0] = t00;
            b_T[(j + 0) * k + i + 1] = t10;
            b_T[(j + 0) * k + i + 2] = t20;
            b_T[(j + 0) * k + i + 3] = t30;
            b_T[(j + 1) * k + i + 0] = t01;
            b_T[(j + 1) * k + i + 1] = t11;
            b_T[(j + 1) * k + i + 2] = t21;
            b_T[(j + 1) * k + i + 3] = t31;
            b_T[(j + 2) * k + i + 0] = t02;
            b_T[(j + 2) * k + i + 1] = t12;
            b_T[(j + 2) * k + i + 2] = t22;
            b_T[(j + 2) * k + i + 3] = t32;
            b_T[(j + 3) * k + i + 0] = t03;
            b_T[(j + 3) * k + i + 1] = t13;
            b_T[(j + 3) * k + i + 2] = t23;
            b_T[(j + 3) * k + i + 3] = t33;
        }
    }
}

// 4x4 output tile. __effcc_ignore_memory_order: A and B_T are read-only,
// accumulators fully independent. Assumes m, n divisible by 4, k divisible
// by 2.
__efficient__ void dmm_helper(const int16_t *a, const int16_t *b_T, int32_t *z,
                              uint32_t m, uint32_t k, uint32_t n)
{
    for (uint32_t tr = 0; tr < m; tr += 4)
    {
        for (uint32_t tc = 0; tc < n; tc += 4)
        {
            int32_t sum00 = 0, sum01 = 0, sum02 = 0, sum03 = 0, sum10 = 0,
                    sum11 = 0, sum12 = 0, sum13 = 0, sum20 = 0, sum21 = 0,
                    sum22 = 0, sum23 = 0, sum30 = 0, sum31 = 0, sum32 = 0,
                    sum33 = 0;

            __effcc_ignore_memory_order for (uint32_t i = 0; i < k; i += 2)
            {
                _v2s16_t a0v = {a[(tr + 0) * k + i], a[(tr + 0) * k + i + 1]};
                _v2s16_t a1v = {a[(tr + 1) * k + i], a[(tr + 1) * k + i + 1]};
                _v2s16_t a2v = {a[(tr + 2) * k + i], a[(tr + 2) * k + i + 1]};
                _v2s16_t a3v = {a[(tr + 3) * k + i], a[(tr + 3) * k + i + 1]};

                _v2s16_t b0v = {b_T[(tc + 0) * k + i],
                                b_T[(tc + 0) * k + i + 1]};
                _v2s16_t b1v = {b_T[(tc + 1) * k + i],
                                b_T[(tc + 1) * k + i + 1]};
                _v2s16_t b2v = {b_T[(tc + 2) * k + i],
                                b_T[(tc + 2) * k + i + 1]};
                _v2s16_t b3v = {b_T[(tc + 3) * k + i],
                                b_T[(tc + 3) * k + i + 1]};

                sum00 += __builtin_effcc_vdot_s16(a0v, b0v);
                sum01 += __builtin_effcc_vdot_s16(a0v, b1v);
                sum02 += __builtin_effcc_vdot_s16(a0v, b2v);
                sum03 += __builtin_effcc_vdot_s16(a0v, b3v);
                sum10 += __builtin_effcc_vdot_s16(a1v, b0v);
                sum11 += __builtin_effcc_vdot_s16(a1v, b1v);
                sum12 += __builtin_effcc_vdot_s16(a1v, b2v);
                sum13 += __builtin_effcc_vdot_s16(a1v, b3v);
                sum20 += __builtin_effcc_vdot_s16(a2v, b0v);
                sum21 += __builtin_effcc_vdot_s16(a2v, b1v);
                sum22 += __builtin_effcc_vdot_s16(a2v, b2v);
                sum23 += __builtin_effcc_vdot_s16(a2v, b3v);
                sum30 += __builtin_effcc_vdot_s16(a3v, b0v);
                sum31 += __builtin_effcc_vdot_s16(a3v, b1v);
                sum32 += __builtin_effcc_vdot_s16(a3v, b2v);
                sum33 += __builtin_effcc_vdot_s16(a3v, b3v);
            }

            z[(tr + 0) * n + tc + 0] = sum00;
            z[(tr + 0) * n + tc + 1] = sum01;
            z[(tr + 0) * n + tc + 2] = sum02;
            z[(tr + 0) * n + tc + 3] = sum03;
            z[(tr + 1) * n + tc + 0] = sum10;
            z[(tr + 1) * n + tc + 1] = sum11;
            z[(tr + 1) * n + tc + 2] = sum12;
            z[(tr + 1) * n + tc + 3] = sum13;
            z[(tr + 2) * n + tc + 0] = sum20;
            z[(tr + 2) * n + tc + 1] = sum21;
            z[(tr + 2) * n + tc + 2] = sum22;
            z[(tr + 2) * n + tc + 3] = sum23;
            z[(tr + 3) * n + tc + 0] = sum30;
            z[(tr + 3) * n + tc + 1] = sum31;
            z[(tr + 3) * n + tc + 2] = sum32;
            z[(tr + 3) * n + tc + 3] = sum33;
        }
    }
}

static int16_t *dmmTemp = NULL;
static uint32_t tempSize = 0;

void dmm(const int16_t *a, // m x k
         const int16_t *b, // k x n
         int32_t *z,       // m x n
         uint32_t m, uint32_t k, uint32_t n)
{
    uint32_t reqSize = k * n * sizeof(int16_t);
    if (dmmTemp == NULL || tempSize < reqSize)
    {
        free(dmmTemp);
        dmmTemp = (int16_t *)malloc(reqSize);
        tempSize = reqSize;
    }
    transpose_b(b, dmmTemp, k, n);
    dmm_helper(a, dmmTemp, z, m, k, n);
}
