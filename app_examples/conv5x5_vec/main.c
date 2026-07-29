#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "conv5x5_vec.h"

#define NUM_ITERATIONS 1

#define N 64
#define M 64

#define CONV5x5_N_PAD (N + 4)
#define CONV5x5_RANDOMIZE_FILTER 1
#define CONV5x5_RANGE 8

const int32_t EXPECTED_CHECKSUM = 1311000512;

data_t filter_global[5 * 5] = {
    DATA_INIT(0), DATA_INIT(-1), DATA_INIT(2), DATA_INIT(-2), DATA_INIT(0),
    DATA_INIT(-1), DATA_INIT(-2), DATA_INIT(3), DATA_INIT(-2), DATA_INIT(-1),
    DATA_INIT(0), DATA_INIT(-3), DATA_INIT(5), DATA_INIT(-3), DATA_INIT(0),
    DATA_INIT(-2), DATA_INIT(-3), DATA_INIT(5), DATA_INIT(2), DATA_INIT(1),
    DATA_INIT(-1), DATA_INIT(0), DATA_INIT(3), DATA_INIT(1), DATA_INIT(-1)};

data_t _Alignas(4) in[CONV5x5_N_PAD * CONV5x5_N_PAD];
data_t _Alignas(4) test_out[N * N];

void print_matrix(const data_t *mat, int n, int stride)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%08x ", mat[i * stride + j]);
        }
        printf("\n");
    }
}

int randPrev = 5;
inline int pseudo_rand()
{
    randPrev = (randPrev * 3) & (CONV5x5_RANGE - 1);
    return randPrev;
}

int main()
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            test_out[i + j * N] = DATA_INIT(0);
            in[i + j * N] = DATA_INIT(pseudo_rand());
        }
    }

    if (CONV5x5_RANDOMIZE_FILTER)
    {
        for (int i = 0;
             i < (int)(sizeof(filter_global) / sizeof(filter_global[0])); i++)
        {
            int f = pseudo_rand() - (CONV5x5_RANGE - 1) / 2;
            filter_global[i] = DATA_INIT(f);
        }
    }

    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
        dconv5x5(in, N, N, filter_global, test_out);

    // Compute scalar reference using strided column-major indexing:
    // out[i + j*n] for i=output_col 0..N-5, j=output_row 0..N-5.
    data_t ref_out[N * N];
    for (int idx = 0; idx < N * N; idx++)
        ref_out[idx] = DATA_INIT(0);
    for (int i = 0; i <= N - 5; i++)
    {
        for (int j = 0; j <= N - 5; j++)
        {
            data_t s = DATA_INIT(0);
            for (int fi = 0; fi < 5; fi++)
                for (int fj = 0; fj < 5; fj++)
                    s += filter_global[fi + fj * 5] *
                         in[(i + fi) + (j + fj) * N];
            ref_out[i + j * N] = s;
        }
    }

    for (int i = 0; i <= N - 5; i++)
    {
        for (int j = 0; j <= N - 5; j++)
        {
            if (test_out[i + j * N] != ref_out[i + j * N])
            {
                printf("[conv5x5] FAIL at col=%d row=%d: got %d expected %d\n",
                       i, j, (int)test_out[i + j * N], (int)ref_out[i + j * N]);
                return 1;
            }
        }
    }

    printf("[conv5x5] PASS\n");
    return 0;
}
