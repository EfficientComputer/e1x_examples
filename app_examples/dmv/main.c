#include "mat.h"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "dmv.h"

#define NUM_ITERATIONS 1

const int32_t EXPECTED_CHECKSUM = 715;

int main()
{
    int32_t z[MAT_REF_SIZE];
    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
        dmv((int32_t *)mat_a, (int32_t *)vector, z, MAT_REF_SIZE, MAT_REF_SIZE);

    int32_t checksum = 0;
    for (int i = 0; i < MAT_REF_SIZE; i++)
    {
        checksum ^= z[i];
    }

    if (checksum == EXPECTED_CHECKSUM)
    {
        printf("[dmv] PASS\n");
    }
    else
    {
        printf("[dmv] FAIL -- checksum %d != %d\n", checksum,
               EXPECTED_CHECKSUM);
    }

    return 0;
}
