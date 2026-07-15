#include "mat.h"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "smv.h"

#define NUM_ITERATIONS 1

#define MAT_REF_SIZE 32
const int32_t EXPECTED_CHECKSUM = 170;

int32_t vector_z[MAT_REF_SIZE];

int main()
{
    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
        smv(mat_a_sparse_indptr, mat_a_sparse_indices, mat_a_sparse_data,
            MAT_REF_SIZE, (const int32_t *)vector, vector_z);

    int32_t checksum = 0;
    for (uint32_t i = 0; i < MAT_REF_SIZE; i++)
    {
        checksum ^= vector_z[i];
    }

    if (checksum == EXPECTED_CHECKSUM)
    {
        printf("[smv] PASS\n");
    }
    else
    {
        printf("[smv] FAIL -- checksum %d != %d\n", checksum,
               EXPECTED_CHECKSUM);
    }

    return 0;
}
