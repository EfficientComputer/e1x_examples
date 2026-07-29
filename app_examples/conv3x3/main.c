#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "conv3x3.h"

#define NUM_ITERATIONS 1

const int32_t EXPECTED_CHECKSUM = -724249960;

data_t filter_global[3 * 3] = {DATA_INIT(3), DATA_INIT(1), DATA_INIT(1),
                               DATA_INIT(1), DATA_INIT(4), DATA_INIT(1),
                               DATA_INIT(2), DATA_INIT(1), DATA_INIT(5)};

// 3 padding rows at the beginning of output
data_t buf[M * N];
data_t test_out[M * N];

int main()
{
    data_t *in = buf;

    // initialize input
    for (int i = 0; i < M * N; i++)
    {
        in[i] = DATA_INIT(2);
    }

    // The first three rows are garbage.

    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
        dconv3x3(in, M, N, INSTRIDE, filter_global, test_out - 3 * OUTSTRIDE);

    int32_t checksum = 0;
    for (int i = 0; i < M - 2; i++)
        for (int j = 0; j < N - 2; j++)
            checksum += (int32_t)test_out[j + i * OUTSTRIDE];

    if (checksum != EXPECTED_CHECKSUM)
    {
        printf("[conv3x3] FAIL checksum: %d != %d\n", checksum,
               EXPECTED_CHECKSUM);
        return 1;
    }

    printf("[conv3x3] PASS\n");
    return 0;
}
