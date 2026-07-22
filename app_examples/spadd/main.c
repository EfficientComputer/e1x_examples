#include "mat.h"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#define NUM_ITERATIONS 1

#define MAT_REF_SIZE 32
const int32_t EXPECTED_CHECKSUM = -30;

void spadd(const int32_t *A_vals, const uint32_t *A_rows,
           const uint32_t *A_coords, const int32_t *B_vals,
           const uint32_t *B_rows, const uint32_t *B_coords, int32_t *Z,
           uint32_t rows);

int32_t mat_z[MAT_REF_SIZE][MAT_REF_SIZE];

int main()
{
    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
        spadd(mat_a_sparse_data, mat_a_sparse_indptr, mat_a_sparse_indices,
              mat_b_sparse_data, mat_b_sparse_indptr, mat_b_sparse_indices,
              (int32_t *)mat_z, MAT_REF_SIZE);

    int32_t checksum = 0;
    for (uint32_t i = 0; i < MAT_REF_SIZE; i++)
    {
        for (uint32_t j = 0; j < MAT_REF_SIZE; j++)
        {
            checksum ^= mat_z[i][j];
        }
    }

    if (checksum == EXPECTED_CHECKSUM)
    {
        printf("[spadd] PASS\n");
    }
    else
    {
        printf("[spadd] FAIL -- checksum %d != %d\n", checksum,
               EXPECTED_CHECKSUM);
    }

    return 0;
}
