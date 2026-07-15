#ifdef __EFFCC__
#include <eff.h>
#endif

#include "conv3x3_vec.h"

// vdot_s16 vectorized 3x3 conv: replace 9 scalar MACs per pixel with
// 4 vdot_s16 pairs + 1 scalar (int32->int16 cast inline).
// Keeps the sliding-window column loop; filter pairs pre-packed before rows.
__efficient__ void dconv3x3_internal(const data_t *in, int m, int n,
                                     int mstride, const data_t *filter,
                                     data_t *out)
{
    data_t f00 = filter[0], f10 = filter[1], f20 = filter[2];
    data_t f01 = filter[3], f11 = filter[4], f21 = filter[5];
    data_t f02 = filter[6], f12 = filter[7], f22 = filter[8];

    // Pre-pack the 9 filter values into 4 pairs (f22 stays scalar)
    _v2s16_t fv0 = {(int16_t)f00, (int16_t)f10};
    _v2s16_t fv1 = {(int16_t)f20, (int16_t)f01};
    _v2s16_t fv2 = {(int16_t)f11, (int16_t)f21};
    _v2s16_t fv3 = {(int16_t)f02, (int16_t)f12};

    for (int i = 0; i <= m - 3; i++)
    {
        for (int j = 0; j <= n - 3; j++)
        {
            _v2s16_t iv0 = {(int16_t)in[(i + 0) * mstride + (j + 0)],
                            (int16_t)in[(i + 0) * mstride + (j + 1)]};
            _v2s16_t iv1 = {(int16_t)in[(i + 0) * mstride + (j + 2)],
                            (int16_t)in[(i + 1) * mstride + (j + 0)]};
            _v2s16_t iv2 = {(int16_t)in[(i + 1) * mstride + (j + 1)],
                            (int16_t)in[(i + 1) * mstride + (j + 2)]};
            _v2s16_t iv3 = {(int16_t)in[(i + 2) * mstride + (j + 0)],
                            (int16_t)in[(i + 2) * mstride + (j + 1)]};

            out[i * n + j] =
                __builtin_effcc_vdot_s16(fv0, iv0) +
                __builtin_effcc_vdot_s16(fv1, iv1) +
                __builtin_effcc_vdot_s16(fv2, iv2) +
                __builtin_effcc_vdot_s16(fv3, iv3) +
                (data_t)(int16_t)f22 *
                    (data_t)(int16_t)in[(i + 2) * mstride + (j + 2)];
        }
    }
}

void dconv3x3(const data_t *in, int m, int n, int mstride, const data_t *filter,
              data_t *out)
{
    // main.c calls with out - 3*OUTSTRIDE; add it back so writes land correctly
    dconv3x3_internal(in, m, n, mstride, filter, out + 3 * n);
}
