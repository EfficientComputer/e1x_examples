#include "conv5x5.h"

#define likely(x) __builtin_expect(!!(x), 1)

#ifdef EFF_BLD_HAND_OPTIMIZED
__efficient__ void dconv5x5_fabric(const data_t *in_row_zero, data_t *out_row_neg5,
                                   int n, int tile_rows_plus_4,
                                   data_t f00, data_t f10, data_t f20, data_t f30, data_t f40,
                                   data_t f01, data_t f11, data_t f21, data_t f31, data_t f41,
                                   data_t f02, data_t f12, data_t f22, data_t f32, data_t f42,
                                   data_t f03, data_t f13, data_t f23, data_t f33, data_t f43,
                                   data_t f04, data_t f14, data_t f24, data_t f34, data_t f44)
{

    for (int col = 0; col <= n - 5; col++)
    {
        data_t in00 = DATA_ZERO, in10 = DATA_ZERO, in20 = DATA_ZERO,
               in30 = DATA_ZERO, in40 = DATA_ZERO;
        data_t in01 = DATA_ZERO, in11 = DATA_ZERO, in21 = DATA_ZERO,
               in31 = DATA_ZERO, in41 = DATA_ZERO;
        data_t in02 = DATA_ZERO, in12 = DATA_ZERO, in22 = DATA_ZERO,
               in32 = DATA_ZERO, in42 = DATA_ZERO;
        data_t in03 = DATA_ZERO, in13 = DATA_ZERO, in23 = DATA_ZERO,
               in33 = DATA_ZERO, in43 = DATA_ZERO;

        const data_t *in = in_row_zero;
        data_t *out = out_row_neg5;

        for (int row = 0; row < tile_rows_plus_4; row++)
        {
            data_t in04 = in[0];
            data_t in14 = in[1];
            data_t in24 = in[2];
            data_t in34 = in[3];
            data_t in44 = in[4];

            in += n;
            out += n;

            data_t val =
                in00 * f00 + in10 * f10 + in20 * f20 + in30 * f30 + in40 * f40;
            val +=
                in01 * f01 + in11 * f11 + in21 * f21 + in31 * f31 + in41 * f41;
            val +=
                in02 * f02 + in12 * f12 + in22 * f22 + in32 * f32 + in42 * f42;
            val +=
                in03 * f03 + in13 * f13 + in23 * f23 + in33 * f33 + in43 * f43;
            val +=
                in04 * f04 + in14 * f14 + in24 * f24 + in34 * f34 + in44 * f44;

            // By only gating only the store (vs the calculation for val)
            // we save a whole bunch of steer gates.
            // It is left to the future to gate the computation of val.
            if (likely(row >= 4))
            {
                *out = val;
            }

            // shift window down
            in00 = in01;
            in10 = in11;
            in20 = in21;
            in30 = in31;
            in40 = in41;

            in01 = in02;
            in11 = in12;
            in21 = in22;
            in31 = in32;
            in41 = in42;

            in02 = in03;
            in12 = in13;
            in22 = in23;
            in32 = in33;
            in42 = in43;

            in03 = in04;
            in13 = in14;
            in23 = in24;
            in33 = in34;
            in43 = in44;
        }

        in_row_zero += 1;
        out_row_neg5 += 1;
    }
}

void dconv5x5_internal(const data_t *in_row_zero, data_t *out_row_neg5,
                       int n, int tile_rows_plus_4,
                       const data_t *pfilter)
{
    data_t f00 = pfilter[0];
    data_t f10 = pfilter[1];
    data_t f20 = pfilter[2];
    data_t f30 = pfilter[3];
    data_t f40 = pfilter[4];
    data_t f01 = pfilter[5];
    data_t f11 = pfilter[6];
    data_t f21 = pfilter[7];
    data_t f31 = pfilter[8];
    data_t f41 = pfilter[9];
    data_t f02 = pfilter[10];
    data_t f12 = pfilter[11];
    data_t f22 = pfilter[12];
    data_t f32 = pfilter[13];
    data_t f42 = pfilter[14];
    data_t f03 = pfilter[15];
    data_t f13 = pfilter[16];
    data_t f23 = pfilter[17];
    data_t f33 = pfilter[18];
    data_t f43 = pfilter[19];
    data_t f04 = pfilter[20];
    data_t f14 = pfilter[21];
    data_t f24 = pfilter[22];
    data_t f34 = pfilter[23];
    data_t f44 = pfilter[24];

    dconv5x5_fabric(in_row_zero, out_row_neg5, n, tile_rows_plus_4,
                    f00, f10, f20, f30, f40,
                    f01, f11, f21, f31, f41,
                    f02, f12, f22, f32, f42,
                    f03, f13, f23, f33, f43,
                    f04, f14, f24, f34, f44);
}
#else

__efficient__ void dconv5x5_internal(const data_t *in, int m, int n,
                                     const data_t *filter, data_t *out)
{
    const int d = 5;
    for (int i = 0; i <= m - d; i++)
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

#endif

void dconv5x5(const data_t *in, int m, int n, const data_t *filter, data_t *out)
{
#ifdef EFF_BLD_HAND_OPTIMIZED
    dconv5x5_internal(in, out - 5 * n, n, n + 4, filter);
#else
    dconv5x5_internal(in, m, n, filter, out);
#endif
}
