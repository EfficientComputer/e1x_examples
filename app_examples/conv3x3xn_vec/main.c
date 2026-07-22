#include <stdio.h>
#include "conv3x3xn_vec.h"

#define NUM_ITERATIONS 1

data_t kernels[NCHANNELS * KERNEL_DIM * KERNEL_DIM];
data_t in[M * INSTRIDE];
// The output stride includes two columns padding columns at
// the beginning so we don't have to gate the output writes
// if FASTPATH is enabled
data_t test_out[M * (N + KERNEL_DIM - 1)];
data_t ref_out[M * N];

void initialize_tensors()
{
    // initialize kernels
    for (int i = 0; i < NCHANNELS * KERNEL_DIM * KERNEL_DIM; ++i)
    {
        kernels[i] = i;
    }

    // If FASTPATH is enabled, the first two columns are invalid because we
    // don't gate output write. If FASTPATH is disabled, you can pass
    // `test_out - KERNEL_DIM + 1` to `dconv_opt` and the output will start at
    // the beginning of the array. For simplicity in these benchmarks, we always
    // pass `test_out` to `dconv_opt` and start reading KERNEL_DIM - 1 cols
    // after the start of the output when testing for correctness.
    data_t *test_out_start = test_out + KERNEL_DIM - 1;

    // initialize inputs and outputs
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            ref_out[N * i + j] = 0;
            test_out_start[OUTSTRIDE * i + j] = 0;
            in[INSTRIDE * i + j] = N * i + j;
        }
    }
}

void dconv_3x3xn_ref(const data_t *inp, int m, int n, int instride,
                     int outstride, const data_t *kern, int nchannels,
                     data_t *out)
{
    for (int i = 0; i < m - KERNEL_DIM + 1; i++)
    {
        for (int j = 0; j < n - KERNEL_DIM + 1; j++)
        {
            data_t sum = 0;
            for (int k = 0; k < nchannels; ++k)
            {
                for (int fi = 0; fi < KERNEL_DIM; fi++)
                {
                    for (int fj = 0; fj < KERNEL_DIM; fj++)
                    {
                        int row = i + fi;
                        int col = j + fj;
                        sum += kern[KERNEL_DIM * KERNEL_DIM * k +
                                    KERNEL_DIM * fi + fj] *
                               inp[col + row * instride];
                    }
                }
            }
            out[i * outstride + j] = sum;
        }
    }
}

int main()
{
    initialize_tensors();

    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
        dconv_3x3xn(in, M, N, INSTRIDE, OUTSTRIDE, kernels, NCHANNELS,
                    test_out);

    dconv_3x3xn_ref(in, M, N, INSTRIDE, N, kernels, NCHANNELS, ref_out);

    // NOTE: do not add optimized_fabric/optimized_sim targets to this app's
    // CMakeLists without also activating the test_out_start offset below.
    // The hand-optimized path writes valid output starting at column
    // KERNEL_DIM-1 (not 0), so the reference comparison would shift by 2
    // and produce a guaranteed FAIL.
    data_t *test_out_start = test_out;

    // check result against reference
    for (int i = 0; i < M - KERNEL_DIM + 1; i++)
    {
        for (int j = 0; j < N - KERNEL_DIM + 1; j++)
        {
            if ((int)ref_out[j + i * N] !=
                (int)test_out_start[j + i * OUTSTRIDE])
            {
                printf(
                    "[conv3x3_multi] FAIL. Mismatch at idx = %d,%d, ref = %d, "
                    "out = "
                    "%d\n",
                    i, j, (int)ref_out[j + i * N],
                    (int)test_out_start[j + i * OUTSTRIDE]);
                return 1;
            }
        }
    }
    printf("[conv3x3xn] PASS\n");
    return 0;
}
