#ifdef __EFFCC__
#include <eff.h>
#endif

#include <stdio.h>
#include <stdint.h>

#include "biquad_filter.h"

#define FIXED_ROUND(x) \
    (((x) + (1 << (FRACBITS - 1))) >> FRACBITS)

#define SAMP_MAX 32767
#define FRACBITS 15

// The hand-optimized kernel packs sample pairs into uint32 accesses, which the
// Fabric maps into its dataflow efficiently. On the RISC-V scalar core those
// 32-bit accesses would be misaligned against the 2-byte-aligned buffers, so
// the scalar build uses the plain per-sample reference kernel below instead.
#if defined(EFF_BLD_HAND_OPTIMIZED) && !defined(EFF_BLD_SCALAR)
__efficient__ void biquad_filter_q15(const int16_t *input, int16_t *output,
                                     int length, const int32_t p1,
                                     const int32_t p2, const int32_t z0,
                                     const int32_t z1, const int32_t z2)
{
    int32_t x1_0 = 0;
    int32_t x2_0 = 0;
    int32_t y1_0 = 0;
    int32_t y2_0 = 0;

    uint32_t in01, in23;
    uint32_t out01, out23;

    int32_t x_0, acc_0, y_0;
    int32_t x2_1, x1_1, y2_1, y1_1;

    int32_t x_1, acc_1, y_1;
    int32_t x2_2, x1_2, y2_2, y1_2;

    int32_t x_2, acc_2, y_2;
    int32_t x2_3, x1_3, y2_3, y1_3;

    int32_t x_3, acc_3, y_3;
    int32_t x2_4, x1_4, y2_4, y1_4;

    const uint32_t *input32 = (const uint32_t *)input;
    uint32_t *output32 = (uint32_t *)output;

    for (int n = 0; n < length; n += 4)
    {
        in01 = input32[0];
        in23 = input32[1];

        x_0 = (int32_t)(int16_t)(in01 & 0xFFFF);
        x_1 = (int32_t)(int16_t)(in01 >> 16);
        x_2 = (int32_t)(int16_t)(in23 & 0xFFFF);
        x_3 = (int32_t)(int16_t)(in23 >> 16);

        // ---- sample n ----
        acc_0 = x_0 * z0 + x1_0 * z1 + x2_0 * z2 - y1_0 * p1 - y2_0 * p2;
        y_0 = FIXED_ROUND(acc_0);

        x2_1 = x1_0;
        x1_1 = x_0;
        y2_1 = y1_0;
        y1_1 = y_0;

        // ---- sample n+1 ----
        acc_1 = x_1 * z0 + x1_1 * z1 + x2_1 * z2 - y1_1 * p1 - y2_1 * p2;
        y_1 = FIXED_ROUND(acc_1);

        x2_2 = x1_1;
        x1_2 = x_1;
        y2_2 = y1_1;
        y1_2 = y_1;

        out01 = ((uint32_t)(uint16_t)y_0) | ((uint32_t)(uint16_t)y_1 << 16);
        output32[0] = out01;

        // ---- sample n+2 ----
        acc_2 = x_2 * z0 + x1_2 * z1 + x2_2 * z2 - y1_2 * p1 - y2_2 * p2;
        y_2 = FIXED_ROUND(acc_2);

        x2_3 = x1_2;
        x1_3 = x_2;
        y2_3 = y1_2;
        y1_3 = y_2;

        // ---- sample n+3 ----
        acc_3 = x_3 * z0 + x1_3 * z1 + x2_3 * z2 - y1_3 * p1 - y2_3 * p2;
        y_3 = FIXED_ROUND(acc_3);

        x2_4 = x1_3;
        x1_4 = x_3;
        y2_4 = y1_3;
        y1_4 = y_3;

        out23 = ((uint32_t)(uint16_t)y_2) | ((uint32_t)(uint16_t)y_3 << 16);
        output32[1] = out23;

        x1_0 = x1_4;
        x2_0 = x2_4;
        y1_0 = y1_4;
        y2_0 = y2_4;

        input32 += 2;
        output32 += 2;
    }
}
#else
__efficient__ void biquad_filter_q15(const int16_t *input, int16_t *output,
                                     int length, const int32_t p1,
                                     const int32_t p2, const int32_t z0,
                                     const int32_t z1, const int32_t z2)
{
    // Implemented from https://en.wikipedia.org/wiki/Digital_biquad_filter

    int32_t x1 = 0;
    int32_t x2 = 0;
    int32_t y1 = 0;
    int32_t y2 = 0;

    for (int n = 0; n < length; n++)
    {
        int32_t x = input[n];
        int32_t acc = x * z0 + x1 * z1 + x2 * z2 - y1 * p1 - y2 * p2;

        int32_t y = FIXED_ROUND(acc);
        output[n] = (int16_t)y;
        y2 = y1;
        y1 = y;
        x2 = x1;
        x1 = x;
    }
}

#endif

void biquad_filter(const int16_t *input, int16_t *output, int length, const float p[2], const float z[3])
{
    biquad_filter_q15(input, output, length,
                      (int32_t)(p[0] * SAMP_MAX),
                      (int32_t)(p[1] * SAMP_MAX),
                      (int32_t)(z[0] * SAMP_MAX),
                      (int32_t)(z[1] * SAMP_MAX),
                      (int32_t)(z[2] * SAMP_MAX));
}
