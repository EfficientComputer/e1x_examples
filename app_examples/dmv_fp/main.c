#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "dmv_fp.h"

#define NUM_ITERATIONS 1

#define M 32
#define MSTRIDE 34
#define N 32

int dmv_float_correctness_test()
{
    // choose values such that all intermediates in the calculation
    // can be precisely represented in fp16.
    float_type a[8 * 8];
    float_type b[8];
    float_type z[8];
    float_type zref[8];

    for (int i = 0; i < 8; ++i)
    {
        b[i] = 2;
        zref[i] = 16;
        for (int j = 0; j < 8; ++j)
        {
            a[i * 8 + j] = 1;
        }
    }

    dmv_fp(a, b, z, 8, 8, 8);

    for (int i = 0; i < 8; ++i)
    {
        uint16_t zi_out;
        memcpy(&zi_out, &z[i], sizeof(uint16_t));
        uint16_t zi_ref;
        memcpy(&zi_ref, &zref[i], sizeof(uint16_t));
        if (zi_out != zi_ref)
        {
            printf("Wrong answer. Mismatch at output element %d\n", i);
            printf("[dmv_fp] FAIL\n");
            return -1;
        }
    }
    printf("[dmv_fp] PASS\n");
    return 0;
}

float_type a[M * MSTRIDE];
float_type b[N];
float_type z[N];

void dmv_float_benchmark()
{
    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
        dmv_fp(a, b, z, M, N, MSTRIDE);
}

int main()
{
    int error;
    // error = dmv_float_correctness_test();
    dmv_float_benchmark();
    return error;
}
