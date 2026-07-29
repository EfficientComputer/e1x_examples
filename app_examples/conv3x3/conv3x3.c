#include "conv3x3.h"

#define likely(x) __builtin_expect(!!(x), 1)

// Inner loop will unroll by 3 when M is divisible by 3
// Sizes and strides must be hardcoded to achieve all compiler opts.
#ifdef EFF_BLD_HAND_OPTIMIZED
__efficient__ void dconv3x3_internal(const data_t *in_row_zero,
                                     data_t *out_row_neg3, int m, int n,
                                     int mstride, data_t f00, data_t f10,
                                     data_t f20, data_t f01, data_t f11,
                                     data_t f21, data_t f02, data_t f12,
                                     data_t f22)
{
    // Need to hardcode this otherwise compiler can't reason about unrolling
    // properly
    for (int col = 0; col < N - 2; col++)
    {
        data_t in00 = DATA_ZERO, in10 = DATA_ZERO, in20 = DATA_ZERO,
               in01 = DATA_ZERO, in11 = DATA_ZERO, in21 = DATA_ZERO;

        const data_t *in = in_row_zero + col;
        data_t *out = out_row_neg3 + col;
        for (int row = 0; row < M; row++)
        {
            data_t in02 = in[0];
            data_t in12 = in[1];
            data_t in22 = in[2];

            // Input loads and output stores will always produce some
            // bank conflicts when the inner loop is unrolled. Hard to get
            // around this unless you have ragged strides for `out`.
            // i.e. strides of N+1, N+1, N+1, N+9, N+1, N+1, N+1, N+9...
            in += INSTRIDE;
            out += OUTSTRIDE;

            // compute convolution output for current 3x3 window
            data_t val = in00 * f00 + in10 * f10 + in20 * f20 +
                         in01 * f01 + in11 * f11 + in21 * f21 +
                         in02 * f02 + in12 * f12 + in22 * f22;

            if (likely(row >= 2))
            {
                *out = val;
            }

            // shift window down
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

__efficient__ void dconv3x3_internal(const data_t *in, int m, int n,
                                     int mstride, const data_t *filter,
                                     data_t *out)
{
    const int d = 3;
    for (int i = 0; i <= m - d; i++)
    {
        for (int j = 0; j <= n - d; j++)
        {
            data_t sum = DATA_ZERO;

            for (int fi = 0; fi < d; fi++)
            {
                for (int fj = 0; fj < d; fj++)
                {
                    int row = i + fi;
                    int col = j + fj;

                    sum += filter[fi * d + fj] * in[col + row * mstride];
                }
            }
            out[i * n + j] = sum;
        }
    }
}

#endif

void dconv3x3(const data_t *in, int m, int n, int mstride, const data_t *filter,
              data_t *out)
{
#ifdef EFF_BLD_HAND_OPTIMIZED
    data_t f00 = filter[0 * 3 + 0], f10 = filter[0 * 3 + 1],
           f20 = filter[0 * 3 + 2], f01 = filter[1 * 3 + 0],
           f11 = filter[1 * 3 + 1], f21 = filter[1 * 3 + 2],
           f02 = filter[2 * 3 + 0], f12 = filter[2 * 3 + 1],
           f22 = filter[2 * 3 + 2];

    dconv3x3_internal(in, out, m, n, mstride, f00, f10, f20, f01, f11, f21, f02,
                      f12, f22);
#else
    dconv3x3_internal(in, m, n, mstride, filter, out + 3 * OUTSTRIDE);
#endif
}
