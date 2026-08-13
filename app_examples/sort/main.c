#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "sort.h"

#define NUM_ITERATIONS 1
#define INPUT_SIZE 512

int32_t temp[INPUT_SIZE];
int32_t originalInput[INPUT_SIZE];
int32_t copiedInput[INPUT_SIZE];

__efficient__ void copyInput(int32_t *restrict dst, int32_t *restrict src,
                             int size)
{
    for (int i = 0; i < size; i++)
    {
        dst[i] = src[i];
    }
}

int main()
{
    srand(10);
    for (int i = 0; i < INPUT_SIZE; i++)
    {
        originalInput[i] = (rand() % 10000) - 5000;
    }

    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
    {
        copyInput(copiedInput, originalInput, INPUT_SIZE);
        merge_sort(copiedInput, temp, INPUT_SIZE);
    }

    for (int i = 1; INPUT_SIZE > i; i++)
    {
        if (copiedInput[i] < copiedInput[i - 1])
        {
            printf("[sort] FAIL index=%d - %d > %d\n", i, copiedInput[i],
                   copiedInput[i - 1]);
            return 0;
        }
    }

    printf("[sort] PASS\n");

    return 0;
}
