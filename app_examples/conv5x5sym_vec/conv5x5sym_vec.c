#ifdef __EFFCC__
#include <eff.h>
#endif

#include "conv5x5sym_vec.h"

// Symmetric 5x5 conv exploiting 4-fold symmetry: f[fi][fj] = f[4-fi][fj]
// = f[fi][4-fj] = f[4-fi][4-fj]. Sum symmetric input pairs first (one add
// per pair), then multiply by the shared coefficient once. Result: 9 unique
// multiplications (4 vdot_s16 + 1 scalar) vs 25 in the non-sym kernel.
//
// Requires input values to fit in int13 so that 4-element sums fit in int16.
// Matches original column-major convention: i=output_col, j=output_row.
__efficient__ void dconv5x5sym_inner(const data_t *in_row_zero,
                                     data_t *out_row_neg5, int n,
                                     int tile_rows_plus_4, data_t f00,
                                     data_t f10, data_t f20, data_t f01,
                                     data_t f11, data_t f21, data_t f02,
                                     data_t f12, data_t f22)
{
    // Pre-pack 9 unique coefficients into 4 vdot pairs + 1 scalar
    _v2s16_t kv0 = {(int16_t)f00, (int16_t)f10};
    _v2s16_t kv1 = {(int16_t)f20, (int16_t)f01};
    _v2s16_t kv2 = {(int16_t)f11, (int16_t)f21};
    _v2s16_t kv3 = {(int16_t)f02, (int16_t)f12};
    int16_t kc = (int16_t)f22;

    for (int i = 0; i <= n - 5; i++)
    {
        for (int j = 0; j + 4 < tile_rows_plus_4; j++)
        {
            // p[fi + fj*n] = in[(i+fi) + (j+fj)*n]  (column-major window)
            const data_t *p = in_row_zero + i + j * n;

            // Sum symmetric input pairs for each of the 9 unique coefficients.
            // f00: corners (fi=0,4 x fj=0,4)
            int16_t s_f00 = (int16_t)((int32_t)p[0 + 0 * n] + p[4 + 0 * n] +
                                      p[0 + 4 * n] + p[4 + 4 * n]);
            // f10: edges (fi=1,3 x fj=0,4)
            int16_t s_f10 = (int16_t)((int32_t)p[1 + 0 * n] + p[3 + 0 * n] +
                                      p[1 + 4 * n] + p[3 + 4 * n]);
            // f20: top/bottom center (fi=2 x fj=0,4)
            int16_t s_f20 = (int16_t)((int32_t)p[2 + 0 * n] + p[2 + 4 * n]);
            // f01: left/right corners (fi=0,4 x fj=1,3)
            int16_t s_f01 = (int16_t)((int32_t)p[0 + 1 * n] + p[4 + 1 * n] +
                                      p[0 + 3 * n] + p[4 + 3 * n]);
            // f11: inner corners (fi=1,3 x fj=1,3)
            int16_t s_f11 = (int16_t)((int32_t)p[1 + 1 * n] + p[3 + 1 * n] +
                                      p[1 + 3 * n] + p[3 + 3 * n]);
            // f21: inner mid edges (fi=2 x fj=1,3)
            int16_t s_f21 = (int16_t)((int32_t)p[2 + 1 * n] + p[2 + 3 * n]);
            // f02: left/right center (fi=0,4 x fj=2)
            int16_t s_f02 = (int16_t)((int32_t)p[0 + 2 * n] + p[4 + 2 * n]);
            // f12: inner center edges (fi=1,3 x fj=2)
            int16_t s_f12 = (int16_t)((int32_t)p[1 + 2 * n] + p[3 + 2 * n]);

            _v2s16_t sv0 = {s_f00, s_f10};
            _v2s16_t sv1 = {s_f20, s_f01};
            _v2s16_t sv2 = {s_f11, s_f21};
            _v2s16_t sv3 = {s_f02, s_f12};

            out_row_neg5[i + (j + 5) * n] =
                __builtin_effcc_vdot_s16(kv0, sv0) +
                __builtin_effcc_vdot_s16(kv1, sv1) +
                __builtin_effcc_vdot_s16(kv2, sv2) +
                __builtin_effcc_vdot_s16(kv3, sv3) +
                (data_t)kc * (data_t)(int16_t)p[2 + 2 * n];
        }
    }
}
