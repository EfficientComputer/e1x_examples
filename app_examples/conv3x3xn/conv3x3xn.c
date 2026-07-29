#include "conv3x3xn.h"

#define likely(x) __builtin_expect(!!(x), 1)

#ifdef EFF_BLD_HAND_OPTIMIZED
__efficient__ void clear_output(data_t *output, int output_size)
{
    for (int i = 0; i < output_size; ++i)
    {
        output[i] = 0;
    }
}
// Hardcoding the sizes and strides of the input allows the compiler
// to unroll more, and avoid software modulo, yielding better performance.
__efficient__ void dconv_3x3xn_internal(const data_t *in, const data_t *kernels,
                                        int nchannels, data_t *out)
{
    for (int k = 0; k < nchannels; ++k)
    {
        const data_t *kernel_start = kernels + KERNEL_DIM * KERNEL_DIM * k;
        data_t k00 = kernel_start[0];
        data_t k01 = kernel_start[1];
        data_t k02 = kernel_start[2];
        data_t k10 = kernel_start[3];
        data_t k11 = kernel_start[4];
        data_t k12 = kernel_start[5];
        data_t k20 = kernel_start[6];
        data_t k21 = kernel_start[7];
        data_t k22 = kernel_start[8];

        data_t in00, in01, in02, in10, in11, in12, in20, in21, in22;

        // row and col loops combined to prevent offloading to scalar core.
        for (int idx = 0; idx < (M - KERNEL_DIM + 1) * N; idx++)
        {
            int i = idx / N;
            int j = idx % N;

            in02 = in[i * INSTRIDE + j];
            in12 = in[(i + 1) * INSTRIDE + j];
            in22 = in[(i + 2) * INSTRIDE + j];

            data_t sum = k00 * in00 + k01 * in01 + k02 * in02 + k10 * in10 +
                         k11 * in11 + k12 * in12 + k20 * in20 + k21 * in21 +
                         k22 * in22;

#if FASTPATH
            out[i * OUTSTRIDE + j] += sum;
#else
            if (likely(j >= 2))
            {
                out[i * OUTSTRIDE + j] += sum;
            }
#endif

            // shift window sideways
            in00 = in01;
            in10 = in11;
            in20 = in21;
            in01 = in02;
            in11 = in12;
            in21 = in22;
        }
    }
}
#else
__efficient__ void dconv_3x3xn_internal(const data_t *in, int m, int n,
                                        int instride, int outstride,
                                        const data_t *kernels, int nchannels,
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
                        sum += kernels[KERNEL_DIM * KERNEL_DIM * k +
                                       KERNEL_DIM * fi + fj] *
                               in[col + row * instride];
                    }
                }
            }
            out[i * outstride + j] = sum;
        }
    }
}
#endif

void dconv_3x3xn(const data_t *in, int m, int n, int instride, int outstride,
                 const data_t *kernels, int nchannels, data_t *out)
{
#ifdef EFF_BLD_HAND_OPTIMIZED
    clear_output(out, m * (n + KERNEL_DIM - 1));
    dconv_3x3xn_internal(in, kernels, nchannels, out);
#else
    dconv_3x3xn_internal(in, m, n, instride, outstride, kernels, nchannels,
                         out);
#endif
}
