#include <stdlib.h>
#include <stdint.h>
#include <eff.h>

#define TILE_SIZE_X 4
#define TILE_SIZE_Y 4

__efficient__ 
void transpose_reformat(const int32_t *restrict input,
                                      int32_t *restrict output, uint32_t rows,
                                      uint32_t cols) {
    for (int y = 0; y < rows; y += 4) {
        __effcc_parallel(2)
        for (int x = 0; x < cols; x++) {
            output[y * cols + x * 4] = input[x * cols + y];
            output[y * cols + x * 4 + 1] = input[x * cols + (y + 1)];
            output[y * cols + x * 4 + 2] = input[x * cols + (y + 2)];
            output[y * cols + x * 4 + 3] = input[x * cols + (y + 3)];
        }
    }
}

__efficient__
void reformat(const int32_t *input, int32_t *output,
                            uint32_t rows, uint32_t cols) {
    for (int y = 0; y < rows; y += 4) {
        __effcc_parallel(2)
        for (int x = 0; x < cols; x++) {
            output[y * cols + x * 4] = input[y * cols + x];
            output[y * cols + x * 4 + 1] = input[(y + 1) * cols + x];
            output[y * cols + x * 4 + 2] = input[(y + 2) * cols + x];
            output[y * cols + x * 4 + 3] = input[(y + 3) * cols + x];
        }
    }
}

__efficient__
void dmm_helper(const int32_t* a,
                int32_t* b,
                int32_t* z,
                uint32_t hiddenLen,
                uint32_t outRows,
                uint32_t outCols) {

    for (int tr = 0; outRows > tr; tr += TILE_SIZE_Y) {
        for (int tc = 0; outCols > tc; tc += TILE_SIZE_X) {

            int32_t sum00 = 0, sum01 = 0, sum02 = 0, sum03 = 0,
                    sum10 = 0, sum11 = 0, sum12 = 0, sum13 = 0,
                    sum20 = 0, sum21 = 0, sum22 = 0, sum23 = 0,
                    sum30 = 0, sum31 = 0, sum32 = 0, sum33 = 0;

            __effcc_parallel(1)
            for (int i = 0; hiddenLen > i; i++) {
                int trHidden = tr * hiddenLen;
                int fourI = i * 4;
                int32_t a0 = a[trHidden + fourI];
                int32_t a1 = a[trHidden + fourI + 1];
                int32_t a2 = a[trHidden + fourI + 2];
                int32_t a3 = a[trHidden + fourI + 3];

                int tcHidden = tc * hiddenLen;
                int32_t b0 = b[tcHidden + fourI];
                int32_t b1 = b[tcHidden + fourI + 1];
                int32_t b2 = b[tcHidden + fourI + 2];
                int32_t b3 = b[tcHidden + fourI + 3];

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

            // saving results
            int tr0OutCols = (tr + 0) * outCols;
            int tr1OutCols = (tr + 1) * outCols;
            int tr2OutCols = (tr + 2) * outCols;
            int tr3OutCols = (tr + 3) * outCols;
            z[tr0OutCols + (tc + 0)] = sum00;
            z[tr0OutCols + (tc + 1)] = sum01;
            z[tr0OutCols + (tc + 2)] = sum02;
            z[tr0OutCols + (tc + 3)] = sum03;
            z[tr1OutCols + (tc + 0)] = sum10;
            z[tr1OutCols + (tc + 1)] = sum11;
            z[tr1OutCols + (tc + 2)] = sum12;
            z[tr1OutCols + (tc + 3)] = sum13;
            z[tr2OutCols + (tc + 0)] = sum20;
            z[tr2OutCols + (tc + 1)] = sum21;
            z[tr2OutCols + (tc + 2)] = sum22;
            z[tr2OutCols + (tc + 3)] = sum23;
            z[tr3OutCols + (tc + 0)] = sum30;
            z[tr3OutCols + (tc + 1)] = sum31;
            z[tr3OutCols + (tc + 2)] = sum32;
            z[tr3OutCols + (tc + 3)] = sum33;
        }
    }
}

int32_t* dmmTemp = NULL;
uint32_t tempSize = 0;

void dmm(const int32_t* a, // m x k
         const int32_t* b, // k x n
         int32_t *z, // m x n
         uint32_t m, uint32_t k, uint32_t n) {

    uint32_t aPrimeSize = m * k;
    uint32_t bTransposeSize = k * n;

    uint32_t reqTempSize = (aPrimeSize + bTransposeSize) * sizeof(int32_t);

    if (dmmTemp == NULL || tempSize < reqTempSize) {
        if (!dmmTemp) {
            free(dmmTemp);
        }

        dmmTemp = (int32_t*)malloc(reqTempSize);
        tempSize = reqTempSize;
    }

    transpose_reformat(b, dmmTemp, k, n);

    int32_t* aPrime = dmmTemp + bTransposeSize;
    reformat(a, aPrime, m, k);

    dmm_helper(aPrime, dmmTemp, z, k, m, n);
}

__efficient__
void dmm_reference(const int32_t *a,  // n x m
                   const int32_t *b,  // m x o
                   int32_t *z,        // n x o
                   uint32_t n, uint32_t m, uint32_t o) {
    for (uint32_t i = 0; i < n; i++) {
        for (uint32_t k = 0; k < o; k++) {
            int32_t sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;
            for (uint32_t j = 0; j < m; j++) {
                sum0 += a[i * n + j] * b[j * m + (k + 0)];
            }

            z[i * n + (k + 0)] = sum0;
        }
    }
}

