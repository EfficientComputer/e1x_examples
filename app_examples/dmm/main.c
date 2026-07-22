#include "mat.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "dmm.h"

#define NUM_ITERATIONS 1

const int32_t EXPECTED_CHECKSUM = -373;

int32_t _Alignas(16) mat_z[MAT_REF_SIZE][MAT_REF_SIZE];

int main()
{
    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
    {
        dmm((int32_t *)mat_a, (int32_t *)mat_b, (int32_t *)mat_z, MAT_REF_SIZE,
            MAT_REF_SIZE, MAT_REF_SIZE);
    }

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
        printf("PASS\n");
    }
    else
    {
        printf("FAIL -- checksum %d != %d\n", checksum, EXPECTED_CHECKSUM);
    }

    return 0;
}
