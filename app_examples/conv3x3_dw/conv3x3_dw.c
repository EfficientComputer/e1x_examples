#include <stdint.h>
#include <stdio.h>

#include "conv3x3_dw.h"

#ifdef EFF_BLD_HAND_OPTIMIZED
// Depthwise 3x3, stride=1, SAME padding (1,1,1,1)
// N = 1 (single image)
// Input : NHWC int8  [H, W, C]
// Weights: [3,3,C,1] int8  (flattened as ((kh*3+kw)*C + c))
// Output: NHWC int32 [H, W, C]

__efficient__ void dconv_3x3_dw_internal(const int8_t *input,
                                         const int8_t *weights, int H, int W,
                                         int C, int32_t *output)
{
    const size_t rowStride = (size_t)W * (size_t)C;

    for (int c = 0; c < C; ++c)
    {
        // Load 3x3 taps once
        const int8_t *kc = weights + (size_t)c;
        const int32_t k00 = kc[0 * C];
        const int32_t k01 = kc[1 * C];
        const int32_t k02 = kc[2 * C];
        const int32_t k10 = kc[3 * C];
        const int32_t k11 = kc[4 * C];
        const int32_t k12 = kc[5 * C];
        const int32_t k20 = kc[6 * C];
        const int32_t k21 = kc[7 * C];
        const int32_t k22 = kc[8 * C];

        for (int i = 0; i < H; ++i)
        {
            const int row0 = i - 1;
            const int row2 = i + 1;

            // Base pointers at column 0 for this channel c
            const int8_t *r0 =
                (row0 >= 0) ? (input + (size_t)row0 * rowStride + (size_t)c)
                            : NULL;
            const int8_t *r1 =
                (input + (size_t)i * rowStride + (size_t)c); // always valid
            const int8_t *r2 =
                (row2 < H) ? (input + (size_t)row2 * rowStride + (size_t)c)
                           : NULL;

            // Output pointer at (i, 0, c)
            int32_t *o = output + (size_t)i * rowStride + (size_t)c;

            // ----- j = 0 (left pad) -----
            int32_t in00 = 0, in10 = 0, in20 = 0; // left column = padding
            int32_t in01 = 0, in11 = 0, in21 = 0; // center column (col 0)
            if (r0)
                in01 = (int32_t)r0[0 * C];
            in11 = (int32_t)r1[0 * C];
            if (r2)
                in21 = (int32_t)r2[0 * C];

            int32_t in02 = 0, in12 = 0, in22 = 0; // right column (col 1)
            if (r0)
                in02 = (int32_t)r0[1 * C];
            in12 = (int32_t)r1[1 * C];
            if (r2)
                in22 = (int32_t)r2[1 * C];

            *o = k01 * in01 + k02 * in02 + k11 * in11 + k12 * in12 +
                 k21 * in21 + k22 * in22;
            o += C; // -> (i,1,c)

            // Prepare per-row pointers at column 1
            const int8_t *p0 = r0 ? (r0 + C) : NULL;
            const int8_t *p1 = r1 + C; // always valid
            const int8_t *p2 = r2 ? (r2 + C) : NULL;

            // ----- j = 1 .. W-2 (fast inner loop) -----
            for (int j = 1; j + 1 < W; ++j)
            {
                // Load ONLY the new rightmost column (j+1)
                const int32_t new02 = r0 ? (int32_t)*(p0 + C) : 0;
                const int32_t new12 = (int32_t)*(p1 + C);
                const int32_t new22 = r2 ? (int32_t)*(p2 + C) : 0;

                // ---- Slide window to column j (systolic shift) ----
                in00 = in01;
                in10 = in11;
                in20 = in21;
                in01 = in02;
                in11 = in12;
                in21 = in22;
                in02 = new02;
                in12 = new12;
                in22 = new22;

                // Now compute at column j
                *o = k00 * in00 + k01 * in01 + k02 * in02 + k10 * in10 +
                     k11 * in11 + k12 * in12 + k20 * in20 + k21 * in21 +
                     k22 * in22;
                o += C;

                // advance per-row pointers by +C (next column)
                if (r0)
                    p0 += C;
                p1 += C;
                if (r2)
                    p2 += C;
            }

            // ----- j = W-1 (right pad) -----
            // one more slide with padding as the new rightmost column
            // however, no need to slide the delay buffer, instead use the old
            // data in the different dot product position since this is the last
            // dot product in this H-loop
            *o = k00 * in01 + k01 * in02 + k10 * in11 + k11 * in12 +
                 k20 * in21 + k21 * in22;
        }
    }
}
#else

// Naive depthwise 3x3, SAME padding (1,1,1,1), stride=1, N=1
// input  : NHWC int8  [H, W, C]
// weights: [3,3,C,1] int8  flattened as ((ky*3 + kx)*C + c)
// output : NHWC int32 [H, W, C]
__efficient__ void dconv_3x3_dw_internal(const int8_t *input,
                                         const int8_t *weights, int H, int W,
                                         int C, int32_t *output)
{
    const int in_row_stride = W * C;   // NHWC row stride
    const int out_row_stride = OW * C; // NHWC row stride for output

    for (int oy = 0; oy < OH; ++oy)
    {
        const int in_y0 = oy - 1; // top of 3x3 window
        for (int ox = 0; ox < OW; ++ox)
        {
            const int in_x0 = ox - 1; // left of 3x3 window
            const int out_base = oy * out_row_stride + ox * C;

            for (int c = 0; c < C; ++c)
            {
                int32_t acc = 0;
                const int8_t *wc =
                    weights + c; // channel pointer; taps spaced by +C

                // Sum over 3x3 kernel with zero-padding via bounds checks
                for (int ky = 0; ky < 3; ++ky)
                {
                    const int iy = in_y0 + ky;
                    if (iy < 0 || iy >= H)
                        continue;
                    const int in_y_base = iy * in_row_stride;

                    for (int kx = 0; kx < 3; ++kx)
                    {
                        const int ix = in_x0 + kx;
                        if (ix < 0 || ix >= W)
                            continue;

                        const int in_idx =
                            in_y_base + ix * C + c; // (iy, ix, c)
                        const int w_idx =
                            (ky * 3 + kx) * C; // ((ky,kx), c) via wc

                        acc += (int32_t)wc[w_idx] * (int32_t)input[in_idx];
                    }
                }

                output[out_base + c] = acc;
            }
        }
    }
}
#endif

void dconv_3x3_dw(const int8_t *in, const int8_t *kernels, int H, int W, int C,
                  int32_t *out)
{
    dconv_3x3_dw_internal(in, kernels, H, W, C, out);
}
