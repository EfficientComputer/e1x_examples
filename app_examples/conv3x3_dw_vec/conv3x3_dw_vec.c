#ifdef __EFFCC__
#include <eff.h>
#endif

#include <stdint.h>

#include "conv3x3_dw_vec.h"

// vdot_s16 vectorized depthwise 3x3 conv (NHWC, SAME padding, stride=1).
// Per channel: pre-pack 9 int8 kernel taps as 4 _v2s16_t pairs + 1 scalar.
// Per pixel: 4 vdot_s16 + 1 scalar vs 9 scalar MACs.
__efficient__ void dconv_3x3_dw_internal(const int8_t *input,
                                         const int8_t *weights, int H, int W,
                                         int C, int32_t *output)
{
    const size_t rowStride = (size_t)W * (size_t)C;

    for (int c = 0; c < C; ++c)
    {
        const int8_t *kc = weights + (size_t)c;

        // Pre-pack 9 int8 taps into 4 vdot pairs + 1 scalar (widening to int16)
        _v2s16_t kv0 = {(int16_t)kc[0 * C], (int16_t)kc[1 * C]};
        _v2s16_t kv1 = {(int16_t)kc[2 * C], (int16_t)kc[3 * C]};
        _v2s16_t kv2 = {(int16_t)kc[4 * C], (int16_t)kc[5 * C]};
        _v2s16_t kv3 = {(int16_t)kc[6 * C], (int16_t)kc[7 * C]};
        int16_t k8 = (int16_t)kc[8 * C];

        for (int i = 0; i < H; ++i)
        {
            const int row0 = i - 1;
            const int row2 = i + 1;

            const int8_t *r0 =
                (row0 >= 0) ? (input + (size_t)row0 * rowStride + (size_t)c)
                            : NULL;
            const int8_t *r1 = input + (size_t)i * rowStride + (size_t)c;
            const int8_t *r2 =
                (row2 < H) ? (input + (size_t)row2 * rowStride + (size_t)c)
                           : NULL;

            int32_t *o = output + (size_t)i * rowStride + (size_t)c;

            // j=0: left pad column is zeros
            int32_t in01 = r0 ? (int32_t)r0[0] : 0;
            int32_t in11 = (int32_t)r1[0];
            int32_t in21 = r2 ? (int32_t)r2[0] : 0;
            int32_t in02 = r0 ? (int32_t)r0[C] : 0;
            int32_t in12 = (int32_t)r1[C];
            int32_t in22 = r2 ? (int32_t)r2[C] : 0;

            // First pixel (j=0): only right half of 3x3 (left column is
            // padding)
            *o = (int32_t)(int16_t)kc[1 * C] * (int16_t)in01 +
                 (int32_t)(int16_t)kc[2 * C] * (int16_t)in02 +
                 (int32_t)(int16_t)kc[4 * C] * (int16_t)in11 +
                 (int32_t)(int16_t)kc[5 * C] * (int16_t)in12 +
                 (int32_t)(int16_t)kc[7 * C] * (int16_t)in21 +
                 (int32_t)(int16_t)kc[8 * C] * (int16_t)in22;
            o += C;

            const int8_t *p0 = r0 ? (r0 + C) : NULL;
            const int8_t *p1 = r1 + C;
            const int8_t *p2 = r2 ? (r2 + C) : NULL;

            // j=1..W-2: full 3x3 window — use vdot
            for (int j = 1; j + 1 < W; ++j)
            {
                int32_t in00 = in01, in10 = in11, in20 = in21;
                in01 = in02;
                in11 = in12;
                in21 = in22;
                in02 = r0 ? (int32_t)*(p0 + C) : 0;
                in12 = (int32_t)*(p1 + C);
                in22 = r2 ? (int32_t)*(p2 + C) : 0;

                _v2s16_t iv0 = {(int16_t)in00, (int16_t)in01};
                _v2s16_t iv1 = {(int16_t)in02, (int16_t)in10};
                _v2s16_t iv2 = {(int16_t)in11, (int16_t)in12};
                _v2s16_t iv3 = {(int16_t)in20, (int16_t)in21};

                *o = __builtin_effcc_vdot_s16(kv0, iv0) +
                     __builtin_effcc_vdot_s16(kv1, iv1) +
                     __builtin_effcc_vdot_s16(kv2, iv2) +
                     __builtin_effcc_vdot_s16(kv3, iv3) +
                     (int32_t)k8 * (int32_t)(int16_t)in22;
                o += C;

                if (p0)
                    p0 += C;
                p1 += C;
                if (p2)
                    p2 += C;
            }

            // j=W-1: right pad column is zeros
            *o = (int32_t)(int16_t)kc[0 * C] * (int16_t)in01 +
                 (int32_t)(int16_t)kc[1 * C] * (int16_t)in02 +
                 (int32_t)(int16_t)kc[3 * C] * (int16_t)in11 +
                 (int32_t)(int16_t)kc[4 * C] * (int16_t)in12 +
                 (int32_t)(int16_t)kc[6 * C] * (int16_t)in21 +
                 (int32_t)(int16_t)kc[7 * C] * (int16_t)in22;
        }
    }
}

void dconv_3x3_dw(const int8_t *in, const int8_t *kernels, int H, int W, int C,
                  int32_t *out)
{
    dconv_3x3_dw_internal(in, kernels, H, W, C, out);
}
