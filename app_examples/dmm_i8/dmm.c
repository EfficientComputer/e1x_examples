#ifdef __EFFCC__
#include <eff.h>
#endif

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#define EXTRACT_i8(v, n) ((int8_t)((v >> (8 * (n))) & 0xff))
#define PACK_i8(i0, i1, i2, i3)                                       \
    (((uint32_t)(uint8_t)i0) << 0) | (((uint32_t)(uint8_t)i1) << 8) | \
        (((uint32_t)(uint8_t)i2) << 16) | (((uint32_t)(uint8_t)i3) << 24);
typedef int8_t v4i8 __attribute__((ext_vector_type(4)));

__attribute__((always_inline)) int8_t rescale(int32_t x,
                                              int32_t quantized_multiplier,
                                              int shift)
{
    const int32_t round = ((int32_t)1) << (shift - 1);
    int32_t result = x * quantized_multiplier + round;
    result >>= shift;
    return result;
}

#ifdef EFF_BLD_HAND_OPTIMIZED

__efficient__ void reformat(const uint32_t *restrict input,
                            uint32_t *restrict output, uint32_t rows,
                            uint32_t cols)
{
    int cols4 = cols / 4;
    for (int y = 0; y < rows; y++)
    {
        __effcc_parallel(2) for (int x = 0; x < cols4; x++)
        {
            output[x * rows + y] = input[y * cols4 + x];
        }
    }
}

__efficient__ void transpose_reformat(const uint32_t *input, uint32_t *output,
                                      uint32_t rows, uint32_t cols)
{
    int cols4 = cols / 4;
    for (int y = 0; y < (int32_t)rows; y += 4)
    {
        __effcc_parallel(2) for (int x = 0; x < (int32_t)cols4; ++x)
        {
            uint32_t in0, in1, in2, in3;
            uint32_t out0, out1, out2, out3;
            int8_t t00, t01, t02, t03, t10, t11, t12, t13, t20, t21, t22, t23,
                t30, t31, t32, t33;

            in0 = input[y * cols4 + x];
            in1 = input[(y + 1) * cols4 + x];
            in2 = input[(y + 2) * cols4 + x];
            in3 = input[(y + 3) * cols4 + x];
            t00 = EXTRACT_i8(in0, 0);
            t01 = EXTRACT_i8(in1, 0);
            t02 = EXTRACT_i8(in2, 0);
            t03 = EXTRACT_i8(in3, 0);
            t10 = EXTRACT_i8(in0, 1);
            t11 = EXTRACT_i8(in1, 1);
            t12 = EXTRACT_i8(in2, 1);
            t13 = EXTRACT_i8(in3, 1);
            t20 = EXTRACT_i8(in0, 2);
            t21 = EXTRACT_i8(in1, 2);
            t22 = EXTRACT_i8(in2, 2);
            t23 = EXTRACT_i8(in3, 2);
            t30 = EXTRACT_i8(in0, 3);
            t31 = EXTRACT_i8(in1, 3);
            t32 = EXTRACT_i8(in2, 3);
            t33 = EXTRACT_i8(in3, 3);

            output[y * cols4 + x * 4 + 0] = PACK_i8(t00, t01, t02, t03);
            output[y * cols4 + x * 4 + 1] = PACK_i8(t10, t11, t12, t13);
            output[y * cols4 + x * 4 + 2] = PACK_i8(t20, t21, t22, t23);
            output[y * cols4 + x * 4 + 3] = PACK_i8(t30, t31, t32, t33);
        }
    }
}

#define TILE_SIZE_X 4
#define TILE_SIZE_Y 4

__effcc_rip_exact void rescale_matrix(const int32_t *a, uint32_t *b,
                                      uint32_t rows, uint32_t cols,
                                      int32_t multiplier, int shift)
{
    int cols4 = cols / 4;
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols4; ++j)
        {
            int32_t i00 = rescale(a[i * cols + j * 4 + 0], multiplier, shift);
            int32_t i01 = rescale(a[i * cols + j * 4 + 1], multiplier, shift);
            int32_t i02 = rescale(a[i * cols + j * 4 + 2], multiplier, shift);
            int32_t i03 = rescale(a[i * cols + j * 4 + 3], multiplier, shift);
            b[i * cols4 + j] = PACK_i8(i00, i01, i02, i03);
        }
    }
}

__effcc_rip_exact void dmm_helper(const uint32_t *a, uint32_t *b, uint32_t *z,
                                  int32_t *zUnscaled, uint32_t hiddenLen,
                                  uint32_t outRows, uint32_t outCols,
                                  int32_t multiplier, int shift)
{
    int outCols4 = outCols / 4;
    int hiddenLen4 = hiddenLen / 4;
    for (int tr = 0; tr < outRows; tr += TILE_SIZE_Y)
    {
        for (int tc = 0; tc < outCols4; ++tc)
        {
            int32_t sum00 = 0, sum01 = 0, sum02 = 0, sum03 = 0, sum10 = 0,
                    sum11 = 0, sum12 = 0, sum13 = 0, sum20 = 0, sum21 = 0,
                    sum22 = 0, sum23 = 0, sum30 = 0, sum31 = 0, sum32 = 0,
                    sum33 = 0;

            __effcc_parallel(1) for (int i = 0; i < hiddenLen; i++)
            {
                int trHidden = tr * hiddenLen4 + i;
                uint32_t in_a = a[trHidden];
                int32_t a0 = EXTRACT_i8(in_a, 0);
                int32_t a1 = EXTRACT_i8(in_a, 1);
                int32_t a2 = EXTRACT_i8(in_a, 2);
                int32_t a3 = EXTRACT_i8(in_a, 3);

                int tcHidden = tc * hiddenLen + i;
                uint32_t in_b = b[tcHidden];
                int32_t b0 = EXTRACT_i8(in_b, 0);
                int32_t b1 = EXTRACT_i8(in_b, 1);
                int32_t b2 = EXTRACT_i8(in_b, 2);
                int32_t b3 = EXTRACT_i8(in_b, 3);

                // systolic array
                sum00 += a0 * b0;
                sum01 += a0 * b1;
                sum02 += a0 * b2;
                sum03 += a0 * b3;
                sum10 += a1 * b0;
                sum11 += a1 * b1;
                sum12 += a1 * b2;
                sum13 += a1 * b3;
                sum20 += a2 * b0;
                sum21 += a2 * b1;
                sum22 += a2 * b2;
                sum23 += a2 * b3;
                sum30 += a3 * b0;
                sum31 += a3 * b1;
                sum32 += a3 * b2;
                sum33 += a3 * b3;
            }

            int tr0OutCols = (tr + 0) * outCols;
            int tr1OutCols = (tr + 1) * outCols;
            int tr2OutCols = (tr + 2) * outCols;
            int tr3OutCols = (tr + 3) * outCols;
            zUnscaled[tr0OutCols + (4 * tc + 0)] = sum00;
            zUnscaled[tr0OutCols + (4 * tc + 1)] = sum01;
            zUnscaled[tr0OutCols + (4 * tc + 2)] = sum02;
            zUnscaled[tr0OutCols + (4 * tc + 3)] = sum03;
            zUnscaled[tr1OutCols + (4 * tc + 0)] = sum10;
            zUnscaled[tr1OutCols + (4 * tc + 1)] = sum11;
            zUnscaled[tr1OutCols + (4 * tc + 2)] = sum12;
            zUnscaled[tr1OutCols + (4 * tc + 3)] = sum13;
            zUnscaled[tr2OutCols + (4 * tc + 0)] = sum20;
            zUnscaled[tr2OutCols + (4 * tc + 1)] = sum21;
            zUnscaled[tr2OutCols + (4 * tc + 2)] = sum22;
            zUnscaled[tr2OutCols + (4 * tc + 3)] = sum23;
            zUnscaled[tr3OutCols + (4 * tc + 0)] = sum30;
            zUnscaled[tr3OutCols + (4 * tc + 1)] = sum31;
            zUnscaled[tr3OutCols + (4 * tc + 2)] = sum32;
            zUnscaled[tr3OutCols + (4 * tc + 3)] = sum33;
        }
    }
}

int8_t *dmmTemp = NULL;
uint32_t tempSize = 0;

void print_matrix(const int8_t *a, uint32_t n, uint32_t m)
{
    for (uint32_t i = 0; i < m; i++)
    {
        for (uint32_t j = 0; j < n; j++)
        {
            printf("%3d ", a[i * n + j]);
        }
        printf("\n");
    }
}
void print_matrix_32(const int32_t *a, uint32_t n, uint32_t m)
{
    for (uint32_t i = 0; i < m; i++)
    {
        for (uint32_t j = 0; j < n; j++)
        {
            printf("%3d ", a[i * n + j]);
        }
        printf("\n");
    }
}

void dmm(const int8_t *a, // m x k
         const int8_t *b, // k x n
         int8_t *z,       // m x n
         int32_t *zUnscaled, uint32_t m, uint32_t k, uint32_t n,
         int32_t multiplier, int shift)
{
    uint32_t aPrimeSize = m * k;
    uint32_t bTransposeSize = k * n;

    uint32_t reqTempSize = (aPrimeSize + bTransposeSize) * sizeof(int8_t);

    typedef uint32_t data_t;

    if (dmmTemp == NULL || tempSize < reqTempSize)
    {
        if (!dmmTemp)
        {
            free(dmmTemp);
        }

        dmmTemp = (int8_t *)malloc(reqTempSize);
        tempSize = reqTempSize;
    }

    reformat((uint32_t *)b, (uint32_t *)dmmTemp, k, n);

    int8_t *aPrime = dmmTemp + bTransposeSize;
    transpose_reformat((data_t *)a, (data_t *)aPrime, m, k);

    dmm_helper((data_t *)aPrime, (data_t *)dmmTemp, (data_t *)z, zUnscaled, k,
               m, n, multiplier, shift);

    rescale_matrix(zUnscaled, (data_t *)z, m, n, multiplier, shift);
}

#else

__efficient__ void dmm(const int8_t *a, // m x k
                       const int8_t *b, // k x n
                       int8_t *z,       // m x n
                       int32_t *zUnscaled, uint32_t n, uint32_t m, uint32_t o,
                       int32_t multiplier, int shift)
{
    for (uint32_t i = 0; i < n; i++)
    {
        for (uint32_t k = 0; k < o; k++)
        {
            int32_t sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;
            for (uint32_t j = 0; j < m; j++)
            {
                sum0 += a[i * n + j] * b[j * m + (k + 0)];
            }

            zUnscaled[i * n + (k + 0)] = sum0;
            z[i * n + k] = rescale(sum0, multiplier, shift);
        }
    }
}

#endif
