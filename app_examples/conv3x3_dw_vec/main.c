#include <stdio.h>
#include <stdlib.h>

#include "conv3x3_dw_vec.h"

#define NUM_ITERATIONS 1

// Depthwise 3x3 weights: [3,3,C,1] flattened as ((ky*3 + kx)*C + c)
int8_t kernels[KERNEL_DIM * KERNEL_DIM * NCHANNELS];

// Input / bias
int8_t in[HEIGHT * WIDTH * NCHANNELS]; // NHWC
int32_t bias[NCHANNELS];               // bias in accumulator domain (unused in this snippet)

// Outputs
int32_t ref_out[OH * OW * NCHANNELS]; // reference (NHWC, i32)
int8_t test_out[OH * OW * NCHANNELS]; // optional quantized output (NHWC, i8)

int32_t test_out_acc[OH * OW * NCHANNELS]; // optimized output (int32)

void initialize_tensors(void)
{
    // Input [M, N, C] (NHWC): increasing pattern
    int8_t cnt = 0;
    for (int y = 0; y < M; ++y)
    {
        for (int x = 0; x < N; ++x)
        {
            for (int c = 0; c < NCHANNELS; ++c)
            {
                size_t idx =
                    ((size_t)y * (size_t)N + (size_t)x) * (size_t)NCHANNELS +
                    (size_t)c;
                in[idx] = (int8_t)(cnt++);
            }
        }
    }

    // Depthwise weights [3,3,C,1] flattened as ((ky*K + kx)*C + c)
    for (int ky = 0; ky < KERNEL_DIM; ++ky)
    {
        for (int kx = 0; kx < KERNEL_DIM; ++kx)
        {
            const int kpos = ky * KERNEL_DIM + kx; // 0..8
            for (int c = 0; c < NCHANNELS; ++c)
            {
                kernels[(size_t)kpos * (size_t)NCHANNELS + (size_t)c] =
                    (int8_t)(kpos * NCHANNELS + c);
            }
        }
    }
}

void dconv_3x3_dw_ref(const int8_t *input, const int8_t *weights, int H, int W,
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

int main()
{
    initialize_tensors();

    // Run optimized systolic kernel (N=1, NHWC, SAME pad, stride=1)
    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
        dconv_3x3_dw(in,      // NHWC int8
                     kernels, // [3,3,C,1] i8, flattened ((ky*3+kx)*C + c)
                     HEIGHT, WIDTH, NCHANNELS,
                     test_out_acc // NHWC int32
        );

    // Run naive reference (SAME pad, arbitrary STRIDE_{W,H})
    dconv_3x3_dw_ref(in,      // NHWC int8
                     kernels, // [3,3,C,1] i8, flattened ((ky*3+kx)*C + c)
                     HEIGHT, WIDTH, NCHANNELS,
                     ref_out // NHWC int32
    );

    // Compare int32 accumulators
    for (int oy = 0; oy < OH; ++oy)
    {
        for (int ox = 0; ox < OW; ++ox)
        {
            for (int oc = 0; oc < NCHANNELS; ++oc)
            {
                const int idx = ((oy * OW) + ox) * NCHANNELS + oc; // NHWC
                if (ref_out[idx] != test_out_acc[idx])
                {
                    printf(
                        "[conv3x3_dw] FAIL at (oy,ox,oc)=(%d,%d,%d), ref=%d, "
                        "opt=%d\n",
                        oy, ox, oc, (int)ref_out[idx], (int)test_out_acc[idx]);
                    return 1;
                }
            }
        }
    }

    printf("[conv3x3_dw] PASS\n");
    return 0;
}
