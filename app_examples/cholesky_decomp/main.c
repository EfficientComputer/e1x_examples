#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "cholesky.h"

#define NUM_ITERATIONS 1

#define INPUT_SIZE 32

const int32_t expected_checksum = -1391;

int a[INPUT_SIZE * INPUT_SIZE];
int chol[INPUT_SIZE * INPUT_SIZE];

int pseudo_rand(int prev) { return (prev * 1000693) & 0x7FFF; }

int main()
{
    // Init data
    uint16_t prev = 277;
    for (int i = 0; i < INPUT_SIZE; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            prev = pseudo_rand(prev);
            int val = FIXED_POINT_CAST_TO((prev % 7 - 4));
            a[i * INPUT_SIZE + j] = val;
            a[j * INPUT_SIZE + i] = val;
        }
    }

    // make sure A is positive semi-definite
    for (int i = 0; i < INPUT_SIZE; ++i)
    {
        int sum = 0;
        for (int j = 0; j < INPUT_SIZE; ++j)
        {
            if (i != j)
            {
                sum += abs(a[i * INPUT_SIZE + j]);
            }
        }
        a[i * INPUT_SIZE + i] = sum + FIXED_POINT_CAST_TO(1);
    }

    // Run Cholesky
    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
        cholesky(a, INPUT_SIZE, chol, 0);

    // Check
    int32_t checksum = 0;
    for (int i = 0; i < INPUT_SIZE; i++)
    {
        for (int j = 0; j < INPUT_SIZE; j++)
        {
            checksum ^= chol[i * INPUT_SIZE + j];
        }
    }

    if (checksum != expected_checksum)
    {
        printf("[cholesky_decomp] FAIL. Got checksum %d, expected %d\n",
               checksum, expected_checksum);
        return 1;
    }

    printf("[cholesky_decomp] PASS\n");

    return 0;
}
