#ifdef __EFFCC__
#include <eff.h>
#endif

#include "conv5x5_vec.h"

// vdot_s16 vectorized 5x5 conv.
// Matches original column-major convention: i=output_col, j=output_row,
// in[col + row*n], filter[fi + fj*5].
// 12 vdot_s16 pairs + 1 scalar vs 25 scalar MACs per pixel.
__efficient__ void dconv5x5_internal(const data_t *in, int m, int n,
                                     const data_t *filter, data_t *out)
{
    _v2s16_t fv00 = {(int16_t)filter[0], (int16_t)filter[1]};
    _v2s16_t fv01 = {(int16_t)filter[2], (int16_t)filter[3]};
    _v2s16_t fv02 = {(int16_t)filter[4], (int16_t)filter[5]};
    _v2s16_t fv03 = {(int16_t)filter[6], (int16_t)filter[7]};
    _v2s16_t fv04 = {(int16_t)filter[8], (int16_t)filter[9]};
    _v2s16_t fv05 = {(int16_t)filter[10], (int16_t)filter[11]};
    _v2s16_t fv06 = {(int16_t)filter[12], (int16_t)filter[13]};
    _v2s16_t fv07 = {(int16_t)filter[14], (int16_t)filter[15]};
    _v2s16_t fv08 = {(int16_t)filter[16], (int16_t)filter[17]};
    _v2s16_t fv09 = {(int16_t)filter[18], (int16_t)filter[19]};
    _v2s16_t fv10 = {(int16_t)filter[20], (int16_t)filter[21]};
    _v2s16_t fv11 = {(int16_t)filter[22], (int16_t)filter[23]};
    int16_t f24 = (int16_t)filter[24];

    for (int i = 0; i <= m - 5; i++)
    { // i = output column
        for (int j = 0; j <= n - 5; j++)
        { // j = output row
            const data_t *p =
                in + i + j * n; // p[fi + fj*n] = in[(i+fi)+(j+fj)*n]

            _v2s16_t iv00 = {(int16_t)p[0 + 0 * n], (int16_t)p[1 + 0 * n]};
            _v2s16_t iv01 = {(int16_t)p[2 + 0 * n], (int16_t)p[3 + 0 * n]};
            _v2s16_t iv02 = {(int16_t)p[4 + 0 * n], (int16_t)p[0 + 1 * n]};
            _v2s16_t iv03 = {(int16_t)p[1 + 1 * n], (int16_t)p[2 + 1 * n]};
            _v2s16_t iv04 = {(int16_t)p[3 + 1 * n], (int16_t)p[4 + 1 * n]};
            _v2s16_t iv05 = {(int16_t)p[0 + 2 * n], (int16_t)p[1 + 2 * n]};
            _v2s16_t iv06 = {(int16_t)p[2 + 2 * n], (int16_t)p[3 + 2 * n]};
            _v2s16_t iv07 = {(int16_t)p[4 + 2 * n], (int16_t)p[0 + 3 * n]};
            _v2s16_t iv08 = {(int16_t)p[1 + 3 * n], (int16_t)p[2 + 3 * n]};
            _v2s16_t iv09 = {(int16_t)p[3 + 3 * n], (int16_t)p[4 + 3 * n]};
            _v2s16_t iv10 = {(int16_t)p[0 + 4 * n], (int16_t)p[1 + 4 * n]};
            _v2s16_t iv11 = {(int16_t)p[2 + 4 * n], (int16_t)p[3 + 4 * n]};

            out[i + j * n] = __builtin_effcc_vdot_s16(fv00, iv00) +
                             __builtin_effcc_vdot_s16(fv01, iv01) +
                             __builtin_effcc_vdot_s16(fv02, iv02) +
                             __builtin_effcc_vdot_s16(fv03, iv03) +
                             __builtin_effcc_vdot_s16(fv04, iv04) +
                             __builtin_effcc_vdot_s16(fv05, iv05) +
                             __builtin_effcc_vdot_s16(fv06, iv06) +
                             __builtin_effcc_vdot_s16(fv07, iv07) +
                             __builtin_effcc_vdot_s16(fv08, iv08) +
                             __builtin_effcc_vdot_s16(fv09, iv09) +
                             __builtin_effcc_vdot_s16(fv10, iv10) +
                             __builtin_effcc_vdot_s16(fv11, iv11) +
                             (data_t)f24 * (data_t)(int16_t)p[4 + 4 * n];
        }
    }
}

void dconv5x5(const data_t *in, int m, int n, const data_t *filter,
              data_t *out)
{
    dconv5x5_internal(in, m, n, filter, out);
}
