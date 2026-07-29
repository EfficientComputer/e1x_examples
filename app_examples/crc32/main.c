#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "crc32.h"

#define NUM_ITERATIONS 1

#define LENGTH 512
unsigned char message[LENGTH];
unsigned int expected = 0x3d0b414b;

int main()
{
    unsigned int res;

    // Initialize message with pseudo-random data
    unsigned int seed = 123456789;
    for (int i = 0; i < LENGTH - 1; ++i)
    {
        seed = 1664525 * seed + 1013904223;
        message[i] = (char)(seed % 255) + 1; // avoid null byte
    }
    message[LENGTH - 1] = 0;

    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
    {
        crc32b(message, &res);
    }

    printf("[CRC] Expected: %x, Actual: %x -- ", expected, res);
    if (res == expected)
    {
        printf("PASS\n");
    }
    else
    {
        printf("FAIL\n");
        return 1;
    }

    return 0;
}
