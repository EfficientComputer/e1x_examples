#include "mat.h"
#include <stdio.h>
#include <stdint.h>

#include "dmm_vectorzied.h"

#define NUM_ITERATIONS 1

const int32_t EXPECTED_CHECKSUM = -373;

int16_t mat_a_i16[MAT_REF_SIZE][MAT_REF_SIZE] __attribute__((aligned(16)));
int16_t mat_b_i16[MAT_REF_SIZE][MAT_REF_SIZE] __attribute__((aligned(16)));
int32_t mat_z[MAT_REF_SIZE][MAT_REF_SIZE] __attribute__((aligned(16)));

int main()
{
    for (uint32_t i = 0; i < MAT_REF_SIZE; i++)
        for (uint32_t j = 0; j < MAT_REF_SIZE; j++)
            mat_a_i16[i][j] = (int16_t)mat_a[i][j];

    for (uint32_t i = 0; i < MAT_REF_SIZE; i++)
        for (uint32_t j = 0; j < MAT_REF_SIZE; j++)
            mat_b_i16[i][j] = (int16_t)mat_b[i][j];

    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
        dmm((int16_t *)mat_a_i16, (int16_t *)mat_b_i16, (int32_t *)mat_z,
            MAT_REF_SIZE, MAT_REF_SIZE, MAT_REF_SIZE);

    // Element-wise reference matmul (i-k-j order, int16 inputs → int32).
    // XOR fold cannot catch transpose/tiling bugs; explicit reference can.
    int32_t ref_z[MAT_REF_SIZE][MAT_REF_SIZE];
    for (uint32_t i = 0; i < MAT_REF_SIZE; i++)
        for (uint32_t j = 0; j < MAT_REF_SIZE; j++)
            ref_z[i][j] = 0;
    for (uint32_t i = 0; i < MAT_REF_SIZE; i++)
        for (uint32_t k = 0; k < MAT_REF_SIZE; k++)
        {
            int32_t a = mat_a_i16[i][k];
            for (uint32_t j = 0; j < MAT_REF_SIZE; j++)
                ref_z[i][j] += a * (int32_t)mat_b_i16[k][j];
        }

    for (uint32_t i = 0; i < MAT_REF_SIZE; i++)
    {
        for (uint32_t j = 0; j < MAT_REF_SIZE; j++)
        {
            if (mat_z[i][j] != ref_z[i][j])
            {
                printf("FAIL -- mismatch at [%u][%u]: got %d expected %d\n", i,
                       j, mat_z[i][j], ref_z[i][j]);
                return 1;
            }
        }
    }

    printf("PASS\n");
    return 0;
}
