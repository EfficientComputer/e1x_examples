#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "conv5x5sym.h"

#define NUM_ITERATIONS 1

#define N 64

#define N_PAD (N + 4)
#define RANDOMIZE_FILTER 1
#define RANGE 8
#define TILE_ROWS N

data_t filter_global[5 * 5] = {
    DATA_INIT(0), DATA_INIT(-1), DATA_INIT(2), DATA_ZERO, DATA_ZERO,
    DATA_INIT(-1), DATA_INIT(-2), DATA_INIT(3), DATA_ZERO, DATA_ZERO,
    DATA_INIT(0), DATA_INIT(-3), DATA_INIT(5), DATA_ZERO, DATA_ZERO,
    DATA_ZERO, DATA_ZERO, DATA_ZERO, DATA_ZERO, DATA_ZERO,
    DATA_ZERO, DATA_ZERO, DATA_ZERO, DATA_ZERO, DATA_ZERO};

data_t in[N_PAD * N_PAD];
data_t ref_out[N_PAD * N_PAD];
data_t test_out[N * N];

void symmetrize_filter()
{
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            int ii = (i > 2) ? 4 - i : i;
            int jj = (j > 2) ? 4 - j : j;
            filter_global[i * 5 + j] = filter_global[ii * 5 + jj];
        }
    }
}

static void dconv5x5sym(const data_t *in, const data_t *filter, data_t *out,
                        int n, int tile_rows)
{
    // see conv3x3

    assert(n % tile_rows == 0);

    data_t f00 = filter[0 * 5 + 0], f10 = filter[0 * 5 + 1],
           f20 = filter[0 * 5 + 2];
    data_t f01 = filter[1 * 5 + 0], f11 = filter[1 * 5 + 1],
           f21 = filter[1 * 5 + 2];
    data_t f02 = filter[2 * 5 + 0], f12 = filter[2 * 5 + 1],
           f22 = filter[2 * 5 + 2];

    for (int row_start = 0; row_start < n; row_start += tile_rows)
    {
        // check conditions needed for loops
        assert(0 <= n - 5);
        // currently assuming that tile_rows divides n

        dconv5x5sym_inner((data_t *)in + row_start * n,
                          out + (row_start - 5) * n, n, tile_rows + 4, f00, f10,
                          f20, f01, f11, f21, f02, f12, f22);
    }
}

void dconv_ref(const data_t *in, int n, const data_t *filter, int d,
               data_t *out)
{
    for (int i = 0; i <= n - d; i++)
    {
        for (int j = 0; j <= n - d; j++)
        {
            data_t sum = DATA_INIT(0);

            for (int fi = 0; fi < d; fi++)
            {
                for (int fj = 0; fj < d; fj++)
                {
                    int col = i + fi;
                    int row = j + fj;

                    sum += filter[fi + fj * d] * in[col + row * n];
                }
            }

            out[i + j * n] = sum;
        }
    }
}

void print_matrix(const data_t *mat, int n, int stride)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%08x ", mat[i * stride + j][0]);
        }
        printf("\n");
    }
}

int main()
{
    srand(42);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            ref_out[i + j * N] = DATA_INIT(0);
            test_out[i + j * N] = DATA_INIT(0);
            in[i + j * N] = DATA_INIT(rand()) % RANGE;
            /* in[i + j * N] = DATA_INIT(i + j * N); */
        }
    }

    if (RANDOMIZE_FILTER)
    {
        for (int i = 0;
             i < (int)(sizeof(filter_global) / sizeof(filter_global[0])); i++)
        {
            int f = rand() % RANGE - (RANGE - 1) / 2;
            filter_global[i] = DATA_INIT(f);
        }
    }

    symmetrize_filter();

    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
        dconv5x5sym(in, filter_global, test_out, N, TILE_ROWS);

    dconv_ref(in, N, filter_global, 5, ref_out);

    for (int i = 0; i < N - 4; i++)
    {
        for (int j = 0; j < N - 4; j++)
        {
            if ((int)ref_out[i + j * N] != (int)test_out[i + j * N])
            {
                printf("[conv5x5] FAIL\n");
                return 1;
            }
        }
    }

    printf("[conv5x5] PASS\n");
    return 0;
}
