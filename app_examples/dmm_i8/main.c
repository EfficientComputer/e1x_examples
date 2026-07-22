#include "mat.h"
#include <stdio.h>
#include <stdlib.h>

#include "dmm.h"

#define NUM_ITERATIONS 1

const uint32_t EXPECTED_CHECKSUM = 901583695;

int8_t mat_a_i8[MAT_REF_SIZE][MAT_REF_SIZE] __attribute__((aligned(16)));
int8_t mat_b_i8[MAT_REF_SIZE][MAT_REF_SIZE] __attribute__((aligned(16)));
int8_t mat_z_i8[MAT_REF_SIZE][MAT_REF_SIZE] __attribute__((aligned(16)));
int32_t mat_z[MAT_REF_SIZE][MAT_REF_SIZE] __attribute__((aligned(16)));

volatile int32_t multiplier = 1;
volatile int shift = 0;

unsigned int crc32b(unsigned char *message, int size)
{
    unsigned int byte, crc, mask;

    crc = 0xFFFFFFFF;
    for (int i = 0; i < size; ++i)
    {
        byte = message[i];
        crc = crc ^ byte;
        for (int j = 0; j < 8; ++j)
        {
            mask = -(crc & 1);
            crc = (crc >> 1) ^ (0xEDB88320 & mask);
        }
    }
    return ~crc;
}

int main()
{
    for (uint32_t i = 0; i < MAT_REF_SIZE; i++)
        for (uint32_t j = 0; j < MAT_REF_SIZE; j++)
            mat_a_i8[i][j] = mat_a[i][j];

    for (uint32_t i = 0; i < MAT_REF_SIZE; i++)
        for (uint32_t j = 0; j < MAT_REF_SIZE; j++)
            mat_b_i8[i][j] = mat_b[i][j];

    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
        dmm((uint32_t *)mat_a_i8, (uint32_t *)mat_b_i8, (uint32_t *)mat_z_i8,
            (int32_t *)mat_z, MAT_REF_SIZE, MAT_REF_SIZE, MAT_REF_SIZE,
            multiplier, shift);

    uint32_t checksum = 0;
    checksum = crc32b((unsigned char *)mat_z_i8,
                      MAT_REF_SIZE * MAT_REF_SIZE * sizeof(int8_t));

    if (checksum == EXPECTED_CHECKSUM)
    {
        printf("PASS\n");
    }
    else
    {
        printf("FAIL -- checksum %u != %u\n", checksum, EXPECTED_CHECKSUM);
    }

    return 0;
}
