// <<AUTOBENCH>> efficient_e1
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include "effintrinsic.h"

typedef uint32_t idx_t;

typedef int8_t vec4xi8_t __attribute__((ext_vector_type(4)));

extern int32_t sdot_vec4xi8(uint32_t a, uint32_t b);
extern void mul_vec4xi8_retptr_4xi32(uint32_t a, uint32_t b, int32_t *c0,
                                     int32_t *c1, int32_t *c2, int32_t *c3);

typedef union {
    vec4xi8_t vi8;
    int8_t ai8[4];
    uint32_t u32;
} word4b_t;

// ======= zero point folding functions ==============
// Filter: [OC, FH, FW, IC]
// Bias  : [OC]
// NOTE: this is identical to conv2d version, can merge.
__efficient__ void conv1d_bias_zp(const int8_t *filter, const int32_t *bias,
                                  idx_t K, idx_t OC, int32_t zp,
                                  int32_t *bias_w_zp) {
    __effcc_ignore_memory_order {
        for (idx_t oc = 0; oc < OC; ++oc) {
            int32_t bias_zp = 0;
            // __effcc_parallel(16)
            for (idx_t k = 0; k < K; k++) {
                bias_zp += (int32_t)filter[oc * K + k];
            }
            // TODO: figure out if sub or add neg is better
            bias_w_zp[oc] = bias[oc] - zp * bias_zp;
        }
    }
}

// Filter: [FH, FW, OC] depthwise (OC==IC)
// Bias  : [OC]
// this basically just operates on transposed filter tensor if IC=1?
// also same as 2d
__efficient__ void depthwise_conv1d_bias_zp(const int8_t *filter,
                                            const int32_t *bias,
                                            idx_t K,  // FH * FW
                                            idx_t OC, int32_t zp,
                                            int32_t *bias_w_zp) {
    __effcc_ignore_memory_order {
        // __effcc_parallel(16)
        for (idx_t oc = 0; oc < OC; oc++) {
            int32_t bias_zp = 0;
            for (idx_t k = 0; k < K; k++) {
                bias_zp += (int32_t)filter[k * OC + oc];
            }
            bias_w_zp[oc] = bias[oc] - zp * bias_zp;
        }
    }
}
// --------------------------------------------

// ============================================================================
// Layer 0: depthwise_conv2d
//   input:  [N][30][1][40] i8
//   weight: [3][1][40][1]  i8   (KH x KW x IC x M)
//   output: [N][28][1][40] i32
// ============================================================================
__effcc_rip_exact void opt_dw_conv2d_layer_0(
    const int8_t *in,      // [IH:30, OC:40]:int8
    const int8_t *filter,  // [FH: 3, OC:40]:int8
    const int32_t *bias,   // [       OC:40]:int32
    int32_t *out           // [OH:28, OC:40]:int32
) {
    const int IH = 30;
    const int OC = 40;
    const int KH = 3;
    const int OH = IH - KH + 1;  // 28

    // tile size: vectorize OCV channels at a time
    const int OCV = 4;
    const int OCR = 1;
    const int OHR = 1;

    __effcc_ignore_memory_order {
        const int8_t *in_oc = in;
        const int8_t *filter_oc = filter;
        const int32_t *bias_oc = bias;
        int32_t *out_oc = out;
        for (int oc = 0; oc < OC / OCV; oc++) {
            uint32_t f_kh0 = *((const uint32_t *)(filter_oc + 0 * OC));
            uint32_t f_kh1 = *((const uint32_t *)(filter_oc + 1 * OC));
            uint32_t f_kh2 = *((const uint32_t *)(filter_oc + 2 * OC));
            int32_t b0 = bias_oc[0], b1 = bias_oc[1], b2 = bias_oc[2],
                    b3 = bias_oc[3];
            const int8_t *in_oc_oh = in_oc;
            int32_t *out_oc_oh = out_oc;
            for (int oh = 0; oh < OH; oh++) {
                int32_t o_oc0_oh0 = b0, o_oc1_oh0 = b1, o_oc2_oh0 = b2,
                        o_oc3_oh0 = b3;
                int32_t mul_oc0, mul_oc1, mul_oc2, mul_oc3;

                uint32_t i_kh0, i_kh1, i_kh2;

                i_kh0 = *((const uint32_t *)(in_oc_oh + 0 * OC));
                mul_vec4xi8_retptr_4xi32(i_kh0, f_kh0, &mul_oc0, &mul_oc1,
                                         &mul_oc2, &mul_oc3);
                o_oc0_oh0 += mul_oc0;
                o_oc1_oh0 += mul_oc1;
                o_oc2_oh0 += mul_oc2;
                o_oc3_oh0 += mul_oc3;

                i_kh1 = *((const uint32_t *)(in_oc_oh + 1 * OC));
                mul_vec4xi8_retptr_4xi32(i_kh1, f_kh1, &mul_oc0, &mul_oc1,
                                         &mul_oc2, &mul_oc3);
                o_oc0_oh0 += mul_oc0;
                o_oc1_oh0 += mul_oc1;
                o_oc2_oh0 += mul_oc2;
                o_oc3_oh0 += mul_oc3;

                i_kh2 = *((const uint32_t *)(in_oc_oh + 2 * OC));
                mul_vec4xi8_retptr_4xi32(i_kh2, f_kh2, &mul_oc0, &mul_oc1,
                                         &mul_oc2, &mul_oc3);
                o_oc0_oh0 += mul_oc0;
                o_oc1_oh0 += mul_oc1;
                o_oc2_oh0 += mul_oc2;
                o_oc3_oh0 += mul_oc3;

                out_oc_oh[0 * OC + 0] = o_oc0_oh0;
                out_oc_oh[0 * OC + 1] = o_oc1_oh0;
                out_oc_oh[0 * OC + 2] = o_oc2_oh0;
                out_oc_oh[0 * OC + 3] = o_oc3_oh0;

                in_oc_oh += OC;
                out_oc_oh += OC;
            }
            in_oc += OCV;
            filter_oc += OCV;
            bias_oc += OCV;
            out_oc += OCV;
        }
    }
}

void dw_conv2d_layer_0(
    const int8_t *in_base, idx_t in_offset, idx_t in_dim_n, idx_t in_dim_h,
    idx_t in_dim_w, idx_t in_dim_c, idx_t in_sn, idx_t in_sh, idx_t in_sw,
    idx_t in_sc, const int8_t *filter_base, idx_t filter_offset,
    // FH=3, FW=1, OC=40, IC=1
    idx_t fd0, idx_t fd1, idx_t fd2, idx_t fd3, idx_t fs0, idx_t fs1, idx_t fs2,
    idx_t fs3, const int32_t *bias_base, idx_t bias_offset, idx_t bias_dim,
    idx_t bias_stride, const int8_t *input_zp_ptr, idx_t input_zp_offset,
    idx_t izp_dim, idx_t izp_stride, const int8_t *weight_zp_ptr,
    idx_t weight_zp_offset, idx_t wzp_dim, idx_t wzp_stride, int32_t *out_base,
    idx_t out_offset, idx_t out_dim_n, idx_t out_dim_h, idx_t out_dim_w,
    idx_t out_dim_c, idx_t out_sn, idx_t out_sh, idx_t out_sw, idx_t out_sc) {
    const int OC = 40;
    const int K = 3;
    int32_t bias_w_zp[OC];
    // assumes offsets are zero.
    depthwise_conv1d_bias_zp(filter_base, bias_base, K, OC, *input_zp_ptr,
                             bias_w_zp);
    // izp folded, not needed since no padding.
    opt_dw_conv2d_layer_0(in_base, filter_base, bias_w_zp, out_base);
}

// ============================================================================
// Layer 1: conv2d (pointwise 1x1)
//   input:  [N][28][1][40]  i8
//   weight: [128][1][1][40] i8   (OC x KH x KW x IC)
//   output: [N][28][1][128] i32
// ============================================================================
__effcc_rip_exact void opt_conv2d_layer_1(
    const int8_t *in,      // [OH:28, IC:40]:int8
    const int8_t *filter,  // [OC:128, IC:40]:int8
    const int32_t *bias,   // [OC:128]:int32
    int32_t *out           // [OH:28, OC:128]:int32
) {
    const int OH = 28;
    const int IC = 40;
    const int OC = 128;

    // vectorize IC by ICV=4 using sdot, unroll OH by OHR=4, OC by OCR=8
    const int ICV = 4;
    // 4x8   = 32 < 36
    const int OHR = 4;
    const int OCR = 8;

    __effcc_ignore_memory_order {
        int32_t *out_oc = out;
        const int32_t *bias_oc = bias;
        const int8_t *filter_oc = filter;
        for (int oc = 0; oc < (OC / OCR); oc++) {
            const int8_t *in_oh = in;
            int32_t *out_oc_oh = out_oc;
            for (int oh = 0; oh < ((OH + OHR - 1) / OHR); oh++) {
                int32_t acc_p0_c0 = bias_oc[0];
                int32_t acc_p0_c1 = bias_oc[1];
                int32_t acc_p0_c2 = bias_oc[2];
                int32_t acc_p0_c3 = bias_oc[3];
                int32_t acc_p0_c4 = bias_oc[4];
                int32_t acc_p0_c5 = bias_oc[5];
                int32_t acc_p0_c6 = bias_oc[6];
                int32_t acc_p0_c7 = bias_oc[7];

                int32_t acc_p1_c0 = acc_p0_c0;
                int32_t acc_p1_c1 = acc_p0_c1;
                int32_t acc_p1_c2 = acc_p0_c2;
                int32_t acc_p1_c3 = acc_p0_c3;
                int32_t acc_p1_c4 = acc_p0_c4;
                int32_t acc_p1_c5 = acc_p0_c5;
                int32_t acc_p1_c6 = acc_p0_c6;
                int32_t acc_p1_c7 = acc_p0_c7;

                int32_t acc_p2_c0 = acc_p0_c0;
                int32_t acc_p2_c1 = acc_p0_c1;
                int32_t acc_p2_c2 = acc_p0_c2;
                int32_t acc_p2_c3 = acc_p0_c3;
                int32_t acc_p2_c4 = acc_p0_c4;
                int32_t acc_p2_c5 = acc_p0_c5;
                int32_t acc_p2_c6 = acc_p0_c6;
                int32_t acc_p2_c7 = acc_p0_c7;

                int32_t acc_p3_c0 = acc_p0_c0;
                int32_t acc_p3_c1 = acc_p0_c1;
                int32_t acc_p3_c2 = acc_p0_c2;
                int32_t acc_p3_c3 = acc_p0_c3;
                int32_t acc_p3_c4 = acc_p0_c4;
                int32_t acc_p3_c5 = acc_p0_c5;
                int32_t acc_p3_c6 = acc_p0_c6;
                int32_t acc_p3_c7 = acc_p0_c7;

                const int8_t *filter_oc_ic = filter_oc;
                const int8_t *in_oh_ic = in_oh;
                for (int ic = 0; ic < (IC / ICV); ic++) {
                    word4b_t in_p0_v, in_p1_v, in_p2_v, in_p3_v;
                    word4b_t filt_c0_v, filt_c1_v, filt_c2_v, filt_c3_v;
                    word4b_t filt_c4_v, filt_c5_v, filt_c6_v, filt_c7_v;

                    in_p0_v.u32 = *((const uint32_t *)(in_oh_ic + 0 * IC));
                    in_p1_v.u32 = *((const uint32_t *)(in_oh_ic + 1 * IC));
                    in_p2_v.u32 = *((const uint32_t *)(in_oh_ic + 2 * IC));
                    in_p3_v.u32 = *((const uint32_t *)(in_oh_ic + 3 * IC));

                    filt_c0_v.u32 =
                        *((const uint32_t *)(filter_oc_ic + 0 * IC));
                    filt_c1_v.u32 =
                        *((const uint32_t *)(filter_oc_ic + 1 * IC));
                    filt_c2_v.u32 =
                        *((const uint32_t *)(filter_oc_ic + 2 * IC));
                    filt_c3_v.u32 =
                        *((const uint32_t *)(filter_oc_ic + 3 * IC));
                    filt_c4_v.u32 =
                        *((const uint32_t *)(filter_oc_ic + 4 * IC));
                    filt_c5_v.u32 =
                        *((const uint32_t *)(filter_oc_ic + 5 * IC));
                    filt_c6_v.u32 =
                        *((const uint32_t *)(filter_oc_ic + 6 * IC));
                    filt_c7_v.u32 =
                        *((const uint32_t *)(filter_oc_ic + 7 * IC));

                    acc_p0_c0 += sdot_vec4xi8(in_p0_v.u32, filt_c0_v.u32);
                    acc_p1_c0 += sdot_vec4xi8(in_p1_v.u32, filt_c0_v.u32);
                    acc_p2_c0 += sdot_vec4xi8(in_p2_v.u32, filt_c0_v.u32);
                    acc_p3_c0 += sdot_vec4xi8(in_p3_v.u32, filt_c0_v.u32);

                    acc_p0_c1 += sdot_vec4xi8(in_p0_v.u32, filt_c1_v.u32);
                    acc_p1_c1 += sdot_vec4xi8(in_p1_v.u32, filt_c1_v.u32);
                    acc_p2_c1 += sdot_vec4xi8(in_p2_v.u32, filt_c1_v.u32);
                    acc_p3_c1 += sdot_vec4xi8(in_p3_v.u32, filt_c1_v.u32);

                    acc_p0_c2 += sdot_vec4xi8(in_p0_v.u32, filt_c2_v.u32);
                    acc_p1_c2 += sdot_vec4xi8(in_p1_v.u32, filt_c2_v.u32);
                    acc_p2_c2 += sdot_vec4xi8(in_p2_v.u32, filt_c2_v.u32);
                    acc_p3_c2 += sdot_vec4xi8(in_p3_v.u32, filt_c2_v.u32);

                    acc_p0_c3 += sdot_vec4xi8(in_p0_v.u32, filt_c3_v.u32);
                    acc_p1_c3 += sdot_vec4xi8(in_p1_v.u32, filt_c3_v.u32);
                    acc_p2_c3 += sdot_vec4xi8(in_p2_v.u32, filt_c3_v.u32);
                    acc_p3_c3 += sdot_vec4xi8(in_p3_v.u32, filt_c3_v.u32);

                    acc_p0_c4 += sdot_vec4xi8(in_p0_v.u32, filt_c4_v.u32);
                    acc_p1_c4 += sdot_vec4xi8(in_p1_v.u32, filt_c4_v.u32);
                    acc_p2_c4 += sdot_vec4xi8(in_p2_v.u32, filt_c4_v.u32);
                    acc_p3_c4 += sdot_vec4xi8(in_p3_v.u32, filt_c4_v.u32);

                    acc_p0_c5 += sdot_vec4xi8(in_p0_v.u32, filt_c5_v.u32);
                    acc_p1_c5 += sdot_vec4xi8(in_p1_v.u32, filt_c5_v.u32);
                    acc_p2_c5 += sdot_vec4xi8(in_p2_v.u32, filt_c5_v.u32);
                    acc_p3_c5 += sdot_vec4xi8(in_p3_v.u32, filt_c5_v.u32);

                    acc_p0_c6 += sdot_vec4xi8(in_p0_v.u32, filt_c6_v.u32);
                    acc_p1_c6 += sdot_vec4xi8(in_p1_v.u32, filt_c6_v.u32);
                    acc_p2_c6 += sdot_vec4xi8(in_p2_v.u32, filt_c6_v.u32);
                    acc_p3_c6 += sdot_vec4xi8(in_p3_v.u32, filt_c6_v.u32);

                    acc_p0_c7 += sdot_vec4xi8(in_p0_v.u32, filt_c7_v.u32);
                    acc_p1_c7 += sdot_vec4xi8(in_p1_v.u32, filt_c7_v.u32);
                    acc_p2_c7 += sdot_vec4xi8(in_p2_v.u32, filt_c7_v.u32);
                    acc_p3_c7 += sdot_vec4xi8(in_p3_v.u32, filt_c7_v.u32);

                    in_oh_ic += ICV;
                    filter_oc_ic += ICV;
                }  // end ic loop

                out_oc_oh[0 * OC + 0] = acc_p0_c0;
                out_oc_oh[0 * OC + 1] = acc_p0_c1;
                out_oc_oh[0 * OC + 2] = acc_p0_c2;
                out_oc_oh[0 * OC + 3] = acc_p0_c3;
                out_oc_oh[0 * OC + 4] = acc_p0_c4;
                out_oc_oh[0 * OC + 5] = acc_p0_c5;
                out_oc_oh[0 * OC + 6] = acc_p0_c6;
                out_oc_oh[0 * OC + 7] = acc_p0_c7;

                out_oc_oh[1 * OC + 0] = acc_p1_c0;
                out_oc_oh[1 * OC + 1] = acc_p1_c1;
                out_oc_oh[1 * OC + 2] = acc_p1_c2;
                out_oc_oh[1 * OC + 3] = acc_p1_c3;
                out_oc_oh[1 * OC + 4] = acc_p1_c4;
                out_oc_oh[1 * OC + 5] = acc_p1_c5;
                out_oc_oh[1 * OC + 6] = acc_p1_c6;
                out_oc_oh[1 * OC + 7] = acc_p1_c7;

                out_oc_oh[2 * OC + 0] = acc_p2_c0;
                out_oc_oh[2 * OC + 1] = acc_p2_c1;
                out_oc_oh[2 * OC + 2] = acc_p2_c2;
                out_oc_oh[2 * OC + 3] = acc_p2_c3;
                out_oc_oh[2 * OC + 4] = acc_p2_c4;
                out_oc_oh[2 * OC + 5] = acc_p2_c5;
                out_oc_oh[2 * OC + 6] = acc_p2_c6;
                out_oc_oh[2 * OC + 7] = acc_p2_c7;

                out_oc_oh[3 * OC + 0] = acc_p3_c0;
                out_oc_oh[3 * OC + 1] = acc_p3_c1;
                out_oc_oh[3 * OC + 2] = acc_p3_c2;
                out_oc_oh[3 * OC + 3] = acc_p3_c3;
                out_oc_oh[3 * OC + 4] = acc_p3_c4;
                out_oc_oh[3 * OC + 5] = acc_p3_c5;
                out_oc_oh[3 * OC + 6] = acc_p3_c6;
                out_oc_oh[3 * OC + 7] = acc_p3_c7;

                out_oc_oh += OHR * OC;
                in_oh += OHR * IC;
            }  // end oh loop
            filter_oc += OCR * IC;
            bias_oc += OCR;
            out_oc += OCR;
        }  // end oc loop
    }  // end __effcc_ignore_memory_order
}

void conv2d_layer_1(
    const int8_t *in_base, idx_t in_offset, idx_t in_dim_n, idx_t in_dim_h,
    idx_t in_dim_w, idx_t in_dim_c, idx_t in_sn, idx_t in_sh, idx_t in_sw,
    idx_t in_sc, const int8_t *filter_base, idx_t filter_offset, idx_t fd0,
    idx_t fd1, idx_t fd2, idx_t fd3, idx_t fs0, idx_t fs1, idx_t fs2, idx_t fs3,
    const int32_t *bias_base, idx_t bias_offset, idx_t bias_dim,
    idx_t bias_stride, const int8_t *input_zp_ptr, idx_t input_zp_offset,
    idx_t izp_dim, idx_t izp_stride, const int8_t *weight_zp_ptr,
    idx_t weight_zp_offset, idx_t wzp_dim, idx_t wzp_stride, int32_t *out_base,
    idx_t out_offset, idx_t out_dim_n, idx_t out_dim_h, idx_t out_dim_w,
    idx_t out_dim_c, idx_t out_sn, idx_t out_sh, idx_t out_sw, idx_t out_sc) {
    const int OC = 128;
    const int K = 40;
    int32_t bias_w_zp[OC];
    conv1d_bias_zp(filter_base, bias_base, K, OC, *input_zp_ptr, bias_w_zp);
    opt_conv2d_layer_1(in_base, filter_base, bias_w_zp, out_base);
}

// ============================================================================
// Layer 2: depthwise_conv2d
//   input:  [N][28][1][128] i8
//   weight: [5][1][128][1]  i8
//   output: [N][24][1][128] i32
// ============================================================================
__effcc_rip_exact void opt_dw_conv2d_layer_2(
    const int8_t *in,      // [IH:28, OC:128]:int8
    const int8_t *filter,  // [FH: 5, OC:128]:int8
    const int32_t *bias,   // [       OC:128]:int32
    int32_t *out           // [OH:24, OC:128]:int32
) {
    const int IH = 28;
    const int OC = 128;
    const int KH = 5;
    const int OH = IH - KH + 1;  // 24

    // tile size: vectorize OCV channels at a time
    const int OCV = 4;
    const int OHR = 1;

    __effcc_ignore_memory_order {
        const int8_t *in_oc = in;
        const int8_t *filter_oc = filter;
        const int32_t *bias_oc = bias;
        int32_t *out_oc = out;
        for (int oc = 0; oc < OC / OCV; oc++) {
            // load all KH filter values once — loop invariant w.r.t. OH
            uint32_t f_kh0 = *((const uint32_t *)(filter_oc + 0 * OC));
            uint32_t f_kh1 = *((const uint32_t *)(filter_oc + 1 * OC));
            uint32_t f_kh2 = *((const uint32_t *)(filter_oc + 2 * OC));
            uint32_t f_kh3 = *((const uint32_t *)(filter_oc + 3 * OC));
            uint32_t f_kh4 = *((const uint32_t *)(filter_oc + 4 * OC));

            int32_t b0 = bias_oc[0], b1 = bias_oc[1], b2 = bias_oc[2],
                    b3 = bias_oc[3];

            const int8_t *in_oc_oh = in_oc;
            int32_t *out_oc_oh = out_oc;
            for (int oh = 0; oh < OH; oh++) {
                int32_t o_oc0_oh0 = b0, o_oc1_oh0 = b1, o_oc2_oh0 = b2,
                        o_oc3_oh0 = b3;
                int32_t mul_oc0, mul_oc1, mul_oc2, mul_oc3;
                uint32_t i_kh0, i_kh1, i_kh2, i_kh3, i_kh4;

                i_kh0 = *((const uint32_t *)(in_oc_oh + 0 * OC));
                mul_vec4xi8_retptr_4xi32(i_kh0, f_kh0, &mul_oc0, &mul_oc1,
                                         &mul_oc2, &mul_oc3);
                o_oc0_oh0 += mul_oc0;
                o_oc1_oh0 += mul_oc1;
                o_oc2_oh0 += mul_oc2;
                o_oc3_oh0 += mul_oc3;

                i_kh1 = *((const uint32_t *)(in_oc_oh + 1 * OC));
                mul_vec4xi8_retptr_4xi32(i_kh1, f_kh1, &mul_oc0, &mul_oc1,
                                         &mul_oc2, &mul_oc3);
                o_oc0_oh0 += mul_oc0;
                o_oc1_oh0 += mul_oc1;
                o_oc2_oh0 += mul_oc2;
                o_oc3_oh0 += mul_oc3;

                i_kh2 = *((const uint32_t *)(in_oc_oh + 2 * OC));
                mul_vec4xi8_retptr_4xi32(i_kh2, f_kh2, &mul_oc0, &mul_oc1,
                                         &mul_oc2, &mul_oc3);
                o_oc0_oh0 += mul_oc0;
                o_oc1_oh0 += mul_oc1;
                o_oc2_oh0 += mul_oc2;
                o_oc3_oh0 += mul_oc3;

                i_kh3 = *((const uint32_t *)(in_oc_oh + 3 * OC));
                mul_vec4xi8_retptr_4xi32(i_kh3, f_kh3, &mul_oc0, &mul_oc1,
                                         &mul_oc2, &mul_oc3);
                o_oc0_oh0 += mul_oc0;
                o_oc1_oh0 += mul_oc1;
                o_oc2_oh0 += mul_oc2;
                o_oc3_oh0 += mul_oc3;

                i_kh4 = *((const uint32_t *)(in_oc_oh + 4 * OC));
                mul_vec4xi8_retptr_4xi32(i_kh4, f_kh4, &mul_oc0, &mul_oc1,
                                         &mul_oc2, &mul_oc3);
                o_oc0_oh0 += mul_oc0;
                o_oc1_oh0 += mul_oc1;
                o_oc2_oh0 += mul_oc2;
                o_oc3_oh0 += mul_oc3;

                out_oc_oh[0 * OC + 0] = o_oc0_oh0;
                out_oc_oh[0 * OC + 1] = o_oc1_oh0;
                out_oc_oh[0 * OC + 2] = o_oc2_oh0;
                out_oc_oh[0 * OC + 3] = o_oc3_oh0;

                in_oc_oh += OHR * OC;
                out_oc_oh += OHR * OC;
            }
            in_oc += OCV;
            filter_oc += OCV;
            bias_oc += OCV;
            out_oc += OCV;
        }
    }
}

void dw_conv2d_layer_2(
    const int8_t *in_base, idx_t in_offset, idx_t in_dim_n, idx_t in_dim_h,
    idx_t in_dim_w, idx_t in_dim_c, idx_t in_sn, idx_t in_sh, idx_t in_sw,
    idx_t in_sc, const int8_t *filter_base, idx_t filter_offset, idx_t fd0,
    idx_t fd1, idx_t fd2, idx_t fd3, idx_t fs0, idx_t fs1, idx_t fs2, idx_t fs3,
    const int32_t *bias_base, idx_t bias_offset, idx_t bias_dim,
    idx_t bias_stride, const int8_t *input_zp_ptr, idx_t input_zp_offset,
    idx_t izp_dim, idx_t izp_stride, const int8_t *weight_zp_ptr,
    idx_t weight_zp_offset, idx_t wzp_dim, idx_t wzp_stride, int32_t *out_base,
    idx_t out_offset, idx_t out_dim_n, idx_t out_dim_h, idx_t out_dim_w,
    idx_t out_dim_c, idx_t out_sn, idx_t out_sh, idx_t out_sw, idx_t out_sc) {
    const int OC = 128;
    const int K = 5;
    int32_t bias_w_zp[OC];
    // assumes offsets are zero.
    depthwise_conv1d_bias_zp(filter_base, bias_base, K, OC, *input_zp_ptr,
                             bias_w_zp);
    opt_dw_conv2d_layer_2(in_base, filter_base, bias_w_zp, out_base);
}

// ============================================================================
// Layer 3: conv2d (pointwise 1x1)
//   input:  [N][24][1][128]  i8
//   weight: [128][1][1][128] i8
//   output: [N][24][1][128]  i32
// ============================================================================
__effcc_rip_exact void opt_conv2d_layer_3(
    const int8_t *in,      // [OH:24, IC:128]:int8
    const int8_t *filter,  // [OC:128, IC:128]:int8
    const int32_t *bias,   // [OC:128]:int32
    int32_t *out           // [OH:24, OC:128]:int32
) {
    const int OH = 24;
    const int IC = 128;
    const int OC = 128;

    // vectorize IC by ICV=4 using sdot, unroll OH by OHR=4, OC by OCR=8
    // OH=24 and OC=128 are evenly divisible by OHR=4 and OCR=8
    const int ICV = 4;
    const int OHR = 4;
    const int OCR = 8;
    const int ICR = 1;

    __effcc_ignore_memory_order {
        int32_t *out_oc = out;
        const int32_t *bias_oc = bias;
        const int8_t *filter_oc = filter;
        for (int oc = 0; oc < (OC / OCR); oc++) {
            const int8_t *in_oh = in;
            int32_t *out_oc_oh = out_oc;
            for (int oh = 0; oh < (OH / OHR); oh++) {
                int32_t acc_p0_c0 = bias_oc[0];
                int32_t acc_p0_c1 = bias_oc[1];
                int32_t acc_p0_c2 = bias_oc[2];
                int32_t acc_p0_c3 = bias_oc[3];
                int32_t acc_p0_c4 = bias_oc[4];
                int32_t acc_p0_c5 = bias_oc[5];
                int32_t acc_p0_c6 = bias_oc[6];
                int32_t acc_p0_c7 = bias_oc[7];

                int32_t acc_p1_c0 = acc_p0_c0;
                int32_t acc_p1_c1 = acc_p0_c1;
                int32_t acc_p1_c2 = acc_p0_c2;
                int32_t acc_p1_c3 = acc_p0_c3;
                int32_t acc_p1_c4 = acc_p0_c4;
                int32_t acc_p1_c5 = acc_p0_c5;
                int32_t acc_p1_c6 = acc_p0_c6;
                int32_t acc_p1_c7 = acc_p0_c7;

                int32_t acc_p2_c0 = acc_p0_c0;
                int32_t acc_p2_c1 = acc_p0_c1;
                int32_t acc_p2_c2 = acc_p0_c2;
                int32_t acc_p2_c3 = acc_p0_c3;
                int32_t acc_p2_c4 = acc_p0_c4;
                int32_t acc_p2_c5 = acc_p0_c5;
                int32_t acc_p2_c6 = acc_p0_c6;
                int32_t acc_p2_c7 = acc_p0_c7;

                int32_t acc_p3_c0 = acc_p0_c0;
                int32_t acc_p3_c1 = acc_p0_c1;
                int32_t acc_p3_c2 = acc_p0_c2;
                int32_t acc_p3_c3 = acc_p0_c3;
                int32_t acc_p3_c4 = acc_p0_c4;
                int32_t acc_p3_c5 = acc_p0_c5;
                int32_t acc_p3_c6 = acc_p0_c6;
                int32_t acc_p3_c7 = acc_p0_c7;

                const int8_t *filter_oc_ic = filter_oc;
                const int8_t *in_oh_ic = in_oh;
                // TODO: try doubling up on loads per loop,
                // then iterating over vals and only doing single math op.
                // does this require that inputs aren't stride 16 banks away?
                // likely, otherwise they would be serialized anyway.
                // but different iterations along ic can happen at the same
                // time, so maybe not?
                // TODO: unroll ICR for loads, add loop for ICR for compute
                for (int ic = 0; ic < (IC / (ICR * ICV)); ic++) {
                    word4b_t in_p0_v0, in_p1_v0, in_p2_v0, in_p3_v0;
                    word4b_t filt_c0_v0, filt_c1_v0, filt_c2_v0, filt_c3_v0;
                    word4b_t filt_c4_v0, filt_c5_v0, filt_c6_v0, filt_c7_v0;

                    in_p0_v0.u32 =
                        *((const uint32_t *)(in_oh_ic + 0 * IC + 0 * ICV));
                    in_p1_v0.u32 =
                        *((const uint32_t *)(in_oh_ic + 1 * IC + 0 * ICV));
                    in_p2_v0.u32 =
                        *((const uint32_t *)(in_oh_ic + 2 * IC + 0 * ICV));
                    in_p3_v0.u32 =
                        *((const uint32_t *)(in_oh_ic + 3 * IC + 0 * ICV));

                    filt_c0_v0.u32 =
                        *((const uint32_t *)(filter_oc_ic + 0 * IC + 0 * ICV));
                    filt_c1_v0.u32 =
                        *((const uint32_t *)(filter_oc_ic + 1 * IC + 0 * ICV));
                    filt_c2_v0.u32 =
                        *((const uint32_t *)(filter_oc_ic + 2 * IC + 0 * ICV));
                    filt_c3_v0.u32 =
                        *((const uint32_t *)(filter_oc_ic + 3 * IC + 0 * ICV));
                    filt_c4_v0.u32 =
                        *((const uint32_t *)(filter_oc_ic + 4 * IC + 0 * ICV));
                    filt_c5_v0.u32 =
                        *((const uint32_t *)(filter_oc_ic + 5 * IC + 0 * ICV));
                    filt_c6_v0.u32 =
                        *((const uint32_t *)(filter_oc_ic + 6 * IC + 0 * ICV));
                    filt_c7_v0.u32 =
                        *((const uint32_t *)(filter_oc_ic + 7 * IC + 0 * ICV));

                    // for (int icr = 0; icr < ICR; icr++)
                    {
                        acc_p0_c0 += sdot_vec4xi8(in_p0_v0.u32, filt_c0_v0.u32);
                        acc_p1_c0 += sdot_vec4xi8(in_p1_v0.u32, filt_c0_v0.u32);
                        acc_p2_c0 += sdot_vec4xi8(in_p2_v0.u32, filt_c0_v0.u32);
                        acc_p3_c0 += sdot_vec4xi8(in_p3_v0.u32, filt_c0_v0.u32);

                        acc_p0_c1 += sdot_vec4xi8(in_p0_v0.u32, filt_c1_v0.u32);
                        acc_p1_c1 += sdot_vec4xi8(in_p1_v0.u32, filt_c1_v0.u32);
                        acc_p2_c1 += sdot_vec4xi8(in_p2_v0.u32, filt_c1_v0.u32);
                        acc_p3_c1 += sdot_vec4xi8(in_p3_v0.u32, filt_c1_v0.u32);

                        acc_p0_c2 += sdot_vec4xi8(in_p0_v0.u32, filt_c2_v0.u32);
                        acc_p1_c2 += sdot_vec4xi8(in_p1_v0.u32, filt_c2_v0.u32);
                        acc_p2_c2 += sdot_vec4xi8(in_p2_v0.u32, filt_c2_v0.u32);
                        acc_p3_c2 += sdot_vec4xi8(in_p3_v0.u32, filt_c2_v0.u32);

                        acc_p0_c3 += sdot_vec4xi8(in_p0_v0.u32, filt_c3_v0.u32);
                        acc_p1_c3 += sdot_vec4xi8(in_p1_v0.u32, filt_c3_v0.u32);
                        acc_p2_c3 += sdot_vec4xi8(in_p2_v0.u32, filt_c3_v0.u32);
                        acc_p3_c3 += sdot_vec4xi8(in_p3_v0.u32, filt_c3_v0.u32);

                        acc_p0_c4 += sdot_vec4xi8(in_p0_v0.u32, filt_c4_v0.u32);
                        acc_p1_c4 += sdot_vec4xi8(in_p1_v0.u32, filt_c4_v0.u32);
                        acc_p2_c4 += sdot_vec4xi8(in_p2_v0.u32, filt_c4_v0.u32);
                        acc_p3_c4 += sdot_vec4xi8(in_p3_v0.u32, filt_c4_v0.u32);

                        acc_p0_c5 += sdot_vec4xi8(in_p0_v0.u32, filt_c5_v0.u32);
                        acc_p1_c5 += sdot_vec4xi8(in_p1_v0.u32, filt_c5_v0.u32);
                        acc_p2_c5 += sdot_vec4xi8(in_p2_v0.u32, filt_c5_v0.u32);
                        acc_p3_c5 += sdot_vec4xi8(in_p3_v0.u32, filt_c5_v0.u32);

                        acc_p0_c6 += sdot_vec4xi8(in_p0_v0.u32, filt_c6_v0.u32);
                        acc_p1_c6 += sdot_vec4xi8(in_p1_v0.u32, filt_c6_v0.u32);
                        acc_p2_c6 += sdot_vec4xi8(in_p2_v0.u32, filt_c6_v0.u32);
                        acc_p3_c6 += sdot_vec4xi8(in_p3_v0.u32, filt_c6_v0.u32);

                        acc_p0_c7 += sdot_vec4xi8(in_p0_v0.u32, filt_c7_v0.u32);
                        acc_p1_c7 += sdot_vec4xi8(in_p1_v0.u32, filt_c7_v0.u32);
                        acc_p2_c7 += sdot_vec4xi8(in_p2_v0.u32, filt_c7_v0.u32);
                        acc_p3_c7 += sdot_vec4xi8(in_p3_v0.u32, filt_c7_v0.u32);
                    }

                    in_oh_ic += ICR * ICV;
                    filter_oc_ic += ICR * ICV;
                    // in_oh_ic     += OHR*ICV;
                    // filter_oc_ic += OCR*ICV;
                }  // end ic loop

                out_oc_oh[0 * OC + 0] = acc_p0_c0;
                out_oc_oh[0 * OC + 1] = acc_p0_c1;
                out_oc_oh[0 * OC + 2] = acc_p0_c2;
                out_oc_oh[0 * OC + 3] = acc_p0_c3;
                out_oc_oh[0 * OC + 4] = acc_p0_c4;
                out_oc_oh[0 * OC + 5] = acc_p0_c5;
                out_oc_oh[0 * OC + 6] = acc_p0_c6;
                out_oc_oh[0 * OC + 7] = acc_p0_c7;

                out_oc_oh[1 * OC + 0] = acc_p1_c0;
                out_oc_oh[1 * OC + 1] = acc_p1_c1;
                out_oc_oh[1 * OC + 2] = acc_p1_c2;
                out_oc_oh[1 * OC + 3] = acc_p1_c3;
                out_oc_oh[1 * OC + 4] = acc_p1_c4;
                out_oc_oh[1 * OC + 5] = acc_p1_c5;
                out_oc_oh[1 * OC + 6] = acc_p1_c6;
                out_oc_oh[1 * OC + 7] = acc_p1_c7;

                out_oc_oh[2 * OC + 0] = acc_p2_c0;
                out_oc_oh[2 * OC + 1] = acc_p2_c1;
                out_oc_oh[2 * OC + 2] = acc_p2_c2;
                out_oc_oh[2 * OC + 3] = acc_p2_c3;
                out_oc_oh[2 * OC + 4] = acc_p2_c4;
                out_oc_oh[2 * OC + 5] = acc_p2_c5;
                out_oc_oh[2 * OC + 6] = acc_p2_c6;
                out_oc_oh[2 * OC + 7] = acc_p2_c7;

                out_oc_oh[3 * OC + 0] = acc_p3_c0;
                out_oc_oh[3 * OC + 1] = acc_p3_c1;
                out_oc_oh[3 * OC + 2] = acc_p3_c2;
                out_oc_oh[3 * OC + 3] = acc_p3_c3;
                out_oc_oh[3 * OC + 4] = acc_p3_c4;
                out_oc_oh[3 * OC + 5] = acc_p3_c5;
                out_oc_oh[3 * OC + 6] = acc_p3_c6;
                out_oc_oh[3 * OC + 7] = acc_p3_c7;

                out_oc_oh += OHR * OC;
                in_oh += OHR * IC;
            }  // end oh loop
            filter_oc += OCR * IC;
            bias_oc += OCR;
            out_oc += OCR;
        }  // end oc loop
    }  // end __effcc_ignore_memory_order
}

void conv2d_layer_3(
    const int8_t *in_base, idx_t in_offset, idx_t in_dim_n, idx_t in_dim_h,
    idx_t in_dim_w, idx_t in_dim_c, idx_t in_sn, idx_t in_sh, idx_t in_sw,
    idx_t in_sc, const int8_t *filter_base, idx_t filter_offset, idx_t fd0,
    idx_t fd1, idx_t fd2, idx_t fd3, idx_t fs0, idx_t fs1, idx_t fs2, idx_t fs3,
    const int32_t *bias_base, idx_t bias_offset, idx_t bias_dim,
    idx_t bias_stride, const int8_t *input_zp_ptr, idx_t input_zp_offset,
    idx_t izp_dim, idx_t izp_stride, const int8_t *weight_zp_ptr,
    idx_t weight_zp_offset, idx_t wzp_dim, idx_t wzp_stride, int32_t *out_base,
    idx_t out_offset, idx_t out_dim_n, idx_t out_dim_h, idx_t out_dim_w,
    idx_t out_dim_c, idx_t out_sn, idx_t out_sh, idx_t out_sw, idx_t out_sc) {
    const int OC = 128;
    const int K = 128;
    int32_t bias_w_zp[OC];
    conv1d_bias_zp(filter_base, bias_base, K, OC, *input_zp_ptr, bias_w_zp);
    opt_conv2d_layer_3(in_base, filter_base, bias_w_zp, out_base);
}

// ============================================================================
// Layer 4: depthwise_conv2d
//   input:  [N][24][1][128] i8
//   weight: [10][1][128][1] i8
//   output: [N][15][1][128] i32
// ============================================================================
__effcc_rip_exact void opt_dw_conv2d_layer_4(
    const int8_t *in,      // [IH:24, OC:128]:int8
    const int8_t *filter,  // [FH:10, OC:128]:int8
    const int32_t *bias,   // [       OC:128]:int32
    int32_t *out           // [OH:15, OC:128]:int32
) {
    const int IH = 24;
    const int OC = 128;
    const int KH = 10;
    const int OH = IH - KH + 1;  // 15

    // tile size: vectorize OCV channels at a time
    const int OCV = 4;
    const int OCR = 2;
    const int OHR = 3;  // OH=15 is divisible by 3
    __effcc_ignore_memory_order {
        const int8_t *in_oc = in;
        const int8_t *filter_oc = filter;
        const int32_t *bias_oc = bias;
        int32_t *out_oc = out;
        for (int oc = 0; oc < OC / (OCR * OCV); oc++) {
            const int8_t *in_oc_oh = in_oc;
            int32_t *out_oc_oh = out_oc;
            for (int oh = 0; oh < OH / OHR; oh++) {
                const int8_t *in_oc_oh_kh = in_oc_oh;
                const int8_t *filter_oc_kh = filter_oc;

                int32_t b0 = bias_oc[0], b1 = bias_oc[1], b2 = bias_oc[2],
                        b3 = bias_oc[3], b4 = bias_oc[4], b5 = bias_oc[5],
                        b6 = bias_oc[6], b7 = bias_oc[7];
                int32_t o_oc0_oh0 = b0, o_oc0_oh1 = b0, o_oc0_oh2 = b0,
                        o_oc1_oh0 = b1, o_oc1_oh1 = b1, o_oc1_oh2 = b1,
                        o_oc2_oh0 = b2, o_oc2_oh1 = b2, o_oc2_oh2 = b2,
                        o_oc3_oh0 = b3, o_oc3_oh1 = b3, o_oc3_oh2 = b3,
                        o_oc4_oh0 = b4, o_oc4_oh1 = b4, o_oc4_oh2 = b4,
                        o_oc5_oh0 = b5, o_oc5_oh1 = b5, o_oc5_oh2 = b5,
                        o_oc6_oh0 = b6, o_oc6_oh1 = b6, o_oc6_oh2 = b6,
                        o_oc7_oh0 = b7, o_oc7_oh1 = b7, o_oc7_oh2 = b7;
                for (int kh = 0; kh < KH; kh++) {
                    // filter shared across all OHR rows
                    // each uint32_t holds 4 int8 filter values for different
                    // channels
                    uint32_t f_oc0, f_oc1;
                    f_oc0 = *((const uint32_t *)(filter_oc_kh + 0 * OCV));
                    f_oc1 = *((const uint32_t *)(filter_oc_kh + 1 * OCV));

                    // input layout is [IH:24, OC:128].
                    // each uint32_t holds 4 int8 input values for different
                    // channels
                    uint32_t i_oh0_oc0, i_oh0_oc1, i_oh1_oc0, i_oh1_oc1,
                        i_oh2_oc0, i_oh2_oc1;
                    i_oh0_oc0 =
                        *((const uint32_t *)(in_oc_oh_kh + 0 * OC + 0 * OCV));
                    i_oh0_oc1 =
                        *((const uint32_t *)(in_oc_oh_kh + 0 * OC + 1 * OCV));

                    i_oh1_oc0 =
                        *((const uint32_t *)(in_oc_oh_kh + 1 * OC + 0 * OCV));
                    i_oh1_oc1 =
                        *((const uint32_t *)(in_oc_oh_kh + 1 * OC + 1 * OCV));

                    i_oh2_oc0 =
                        *((const uint32_t *)(in_oc_oh_kh + 2 * OC + 0 * OCV));
                    i_oh2_oc1 =
                        *((const uint32_t *)(in_oc_oh_kh + 2 * OC + 1 * OCV));

                    int32_t mul_oc0, mul_oc1, mul_oc2, mul_oc3, mul_oc4,
                        mul_oc5, mul_oc6, mul_oc7;

                    // row oh+0
                    mul_vec4xi8_retptr_4xi32(i_oh0_oc0, f_oc0, &mul_oc0,
                                             &mul_oc1, &mul_oc2, &mul_oc3);
                    mul_vec4xi8_retptr_4xi32(i_oh0_oc1, f_oc1, &mul_oc4,
                                             &mul_oc5, &mul_oc6, &mul_oc7);
                    o_oc0_oh0 += mul_oc0;
                    o_oc1_oh0 += mul_oc1;
                    o_oc2_oh0 += mul_oc2;
                    o_oc3_oh0 += mul_oc3;
                    o_oc4_oh0 += mul_oc4;
                    o_oc5_oh0 += mul_oc5;
                    o_oc6_oh0 += mul_oc6;
                    o_oc7_oh0 += mul_oc7;

                    // row oh+1
                    mul_vec4xi8_retptr_4xi32(i_oh1_oc0, f_oc0, &mul_oc0,
                                             &mul_oc1, &mul_oc2, &mul_oc3);
                    mul_vec4xi8_retptr_4xi32(i_oh1_oc1, f_oc1, &mul_oc4,
                                             &mul_oc5, &mul_oc6, &mul_oc7);
                    o_oc0_oh1 += mul_oc0;
                    o_oc1_oh1 += mul_oc1;
                    o_oc2_oh1 += mul_oc2;
                    o_oc3_oh1 += mul_oc3;
                    o_oc4_oh1 += mul_oc4;
                    o_oc5_oh1 += mul_oc5;
                    o_oc6_oh1 += mul_oc6;
                    o_oc7_oh1 += mul_oc7;

                    // row oh+2
                    mul_vec4xi8_retptr_4xi32(i_oh2_oc0, f_oc0, &mul_oc0,
                                             &mul_oc1, &mul_oc2, &mul_oc3);
                    mul_vec4xi8_retptr_4xi32(i_oh2_oc1, f_oc1, &mul_oc4,
                                             &mul_oc5, &mul_oc6, &mul_oc7);
                    o_oc0_oh2 += mul_oc0;
                    o_oc1_oh2 += mul_oc1;
                    o_oc2_oh2 += mul_oc2;
                    o_oc3_oh2 += mul_oc3;
                    o_oc4_oh2 += mul_oc4;
                    o_oc5_oh2 += mul_oc5;
                    o_oc6_oh2 += mul_oc6;
                    o_oc7_oh2 += mul_oc7;

                    in_oc_oh_kh += 1 * OC;
                    filter_oc_kh += 1 * OC;
                }
                // store row oh+0
                out_oc_oh[0 * OC + 0] = o_oc0_oh0;
                out_oc_oh[0 * OC + 1] = o_oc1_oh0;
                out_oc_oh[0 * OC + 2] = o_oc2_oh0;
                out_oc_oh[0 * OC + 3] = o_oc3_oh0;
                out_oc_oh[0 * OC + 4] = o_oc4_oh0;
                out_oc_oh[0 * OC + 5] = o_oc5_oh0;
                out_oc_oh[0 * OC + 6] = o_oc6_oh0;
                out_oc_oh[0 * OC + 7] = o_oc7_oh0;

                // store row oh+1
                out_oc_oh[1 * OC + 0] = o_oc0_oh1;
                out_oc_oh[1 * OC + 1] = o_oc1_oh1;
                out_oc_oh[1 * OC + 2] = o_oc2_oh1;
                out_oc_oh[1 * OC + 3] = o_oc3_oh1;
                out_oc_oh[1 * OC + 4] = o_oc4_oh1;
                out_oc_oh[1 * OC + 5] = o_oc5_oh1;
                out_oc_oh[1 * OC + 6] = o_oc6_oh1;
                out_oc_oh[1 * OC + 7] = o_oc7_oh1;

                // store row oh+2
                out_oc_oh[2 * OC + 0] = o_oc0_oh2;
                out_oc_oh[2 * OC + 1] = o_oc1_oh2;
                out_oc_oh[2 * OC + 2] = o_oc2_oh2;
                out_oc_oh[2 * OC + 3] = o_oc3_oh2;
                out_oc_oh[2 * OC + 4] = o_oc4_oh2;
                out_oc_oh[2 * OC + 5] = o_oc5_oh2;
                out_oc_oh[2 * OC + 6] = o_oc6_oh2;
                out_oc_oh[2 * OC + 7] = o_oc7_oh2;

                in_oc_oh += OHR * OC;
                out_oc_oh += OHR * OC;
            }
            in_oc += OCR * OCV;
            filter_oc += OCR * OCV;
            bias_oc += OCR * OCV;
            out_oc += OCR * OCV;
        }
    }
}

void dw_conv2d_layer_4(
    const int8_t *in_base, idx_t in_offset, idx_t in_dim_n, idx_t in_dim_h,
    idx_t in_dim_w, idx_t in_dim_c, idx_t in_sn, idx_t in_sh, idx_t in_sw,
    idx_t in_sc, const int8_t *filter_base, idx_t filter_offset, idx_t fd0,
    idx_t fd1, idx_t fd2, idx_t fd3, idx_t fs0, idx_t fs1, idx_t fs2, idx_t fs3,
    const int32_t *bias_base, idx_t bias_offset, idx_t bias_dim,
    idx_t bias_stride, const int8_t *input_zp_ptr, idx_t input_zp_offset,
    idx_t izp_dim, idx_t izp_stride, const int8_t *weight_zp_ptr,
    idx_t weight_zp_offset, idx_t wzp_dim, idx_t wzp_stride, int32_t *out_base,
    idx_t out_offset, idx_t out_dim_n, idx_t out_dim_h, idx_t out_dim_w,
    idx_t out_dim_c, idx_t out_sn, idx_t out_sh, idx_t out_sw, idx_t out_sc) {
    const int OC = 128;
    const int K = 10;
    int32_t bias_w_zp[OC];
    depthwise_conv1d_bias_zp(filter_base, bias_base, K, OC, *input_zp_ptr,
                             bias_w_zp);
    opt_dw_conv2d_layer_4(in_base, filter_base, bias_w_zp, out_base);
}

// ============================================================================
// Layer 5: conv2d (pointwise 1x1)
//   input:  [N][15][1][128]  i8
//   weight: [128][1][1][128] i8
//   output: [N][15][1][128]  i32
// ============================================================================
// PASS
// 11851 core cyles, 3907 fabric cycles
//  core cycles,  fabric cycles
// instead of 5x7x4, implement 3x6x8; didn't fit
// instead of 5x7x4, implement 3x4x8; 19167         6395
// 3x2x16: 12031         4018
__effcc_rip_exact void opt_conv2d_layer_5(
    const int8_t *in,      // [OH:15, IC:128]:int8
    const int8_t *filter,  // [OC:128, IC:128]:int8
    const int32_t *bias,   // [OC:128]:int32
    int32_t *out           // [OH:15, OC:128]:int32
) {
    const int OH = 15;
    const int IC = 128;
    const int OC = 128;

    // vectorize IC by ICV=4 using sdot, unroll OH by OHWR=5, OC by OCR=7
    const int ICV = 4;
    const int ICR = 4;
    const int OHWR = 3;
    const int OCR = 2;  // should this be increasedd up to 3?

    __effcc_ignore_memory_order {
        int32_t *out_oc = out;
        const int32_t *bias_oc = bias;
        const int8_t *filter_oc = filter;
        for (int oc = 0; oc < (OC / OCR); oc++) {
            const int8_t *in_oh = in;
            int32_t *out_oc_oh = out_oc;

            for (int oh = 0; oh < (OH / OHWR); oh++) {
                int32_t acc_p0_c0 = bias_oc[0];
                int32_t acc_p1_c0 = acc_p0_c0;
                int32_t acc_p2_c0 = acc_p0_c0;
                int32_t acc_p3_c0 = acc_p0_c0;
                int32_t acc_p4_c0 = acc_p0_c0;

                int32_t acc_p0_c1 = bias_oc[1];
                int32_t acc_p1_c1 = acc_p0_c1;
                int32_t acc_p2_c1 = acc_p0_c1;
                int32_t acc_p3_c1 = acc_p0_c1;
                int32_t acc_p4_c1 = acc_p0_c1;

                const int8_t *filter_oc_ic = filter_oc;
                const int8_t *in_oh_ic = in_oh;
                for (int ic = 0; ic < (IC / (ICR * ICV)); ic++) {
                    uint32_t i_oh0_ic0, i_oh1_ic0, i_oh2_ic0, i_oh0_ic1,
                        i_oh1_ic1, i_oh2_ic1, i_oh0_ic2, i_oh1_ic2, i_oh2_ic2,
                        i_oh0_ic3, i_oh1_ic3, i_oh2_ic3;
                    uint32_t f_oc0_ic0, f_oc1_ic0, f_oc2_ic0, f_oc3_ic0,
                        f_oc4_ic0, f_oc5_ic0, f_oc0_ic1, f_oc1_ic1, f_oc2_ic1,
                        f_oc3_ic1, f_oc4_ic1, f_oc5_ic1, f_oc0_ic2, f_oc1_ic2,
                        f_oc2_ic2, f_oc3_ic2, f_oc4_ic2, f_oc5_ic2, f_oc0_ic3,
                        f_oc1_ic3, f_oc2_ic3, f_oc3_ic3, f_oc4_ic3, f_oc5_ic3;

                    i_oh0_ic0 =
                        *((const uint32_t *)(in_oh_ic + 0 * IC + 0 * ICV));
                    i_oh0_ic1 =
                        *((const uint32_t *)(in_oh_ic + 0 * IC + 1 * ICV));
                    i_oh0_ic2 =
                        *((const uint32_t *)(in_oh_ic + 0 * IC + 2 * ICV));
                    i_oh0_ic3 =
                        *((const uint32_t *)(in_oh_ic + 0 * IC + 3 * ICV));

                    i_oh1_ic0 =
                        *((const uint32_t *)(in_oh_ic + 1 * IC + 0 * ICV));
                    i_oh1_ic1 =
                        *((const uint32_t *)(in_oh_ic + 1 * IC + 1 * ICV));
                    i_oh1_ic2 =
                        *((const uint32_t *)(in_oh_ic + 1 * IC + 2 * ICV));
                    i_oh1_ic3 =
                        *((const uint32_t *)(in_oh_ic + 1 * IC + 3 * ICV));

                    i_oh2_ic0 =
                        *((const uint32_t *)(in_oh_ic + 2 * IC + 0 * ICV));
                    i_oh2_ic1 =
                        *((const uint32_t *)(in_oh_ic + 2 * IC + 1 * ICV));
                    i_oh2_ic2 =
                        *((const uint32_t *)(in_oh_ic + 2 * IC + 2 * ICV));
                    i_oh2_ic3 =
                        *((const uint32_t *)(in_oh_ic + 2 * IC + 3 * ICV));

                    f_oc0_ic0 =
                        *((const uint32_t *)(filter_oc_ic + 0 * IC + 0 * ICV));
                    f_oc0_ic1 =
                        *((const uint32_t *)(filter_oc_ic + 0 * IC + 1 * ICV));
                    f_oc0_ic2 =
                        *((const uint32_t *)(filter_oc_ic + 0 * IC + 2 * ICV));
                    f_oc0_ic3 =
                        *((const uint32_t *)(filter_oc_ic + 0 * IC + 3 * ICV));

                    f_oc1_ic0 =
                        *((const uint32_t *)(filter_oc_ic + 1 * IC + 0 * ICV));
                    f_oc1_ic1 =
                        *((const uint32_t *)(filter_oc_ic + 1 * IC + 1 * ICV));
                    f_oc1_ic2 =
                        *((const uint32_t *)(filter_oc_ic + 1 * IC + 2 * ICV));
                    f_oc1_ic3 =
                        *((const uint32_t *)(filter_oc_ic + 1 * IC + 3 * ICV));

                    acc_p0_c0 += sdot_vec4xi8(i_oh0_ic0, f_oc0_ic0);
                    acc_p0_c0 += sdot_vec4xi8(i_oh0_ic1, f_oc0_ic1);
                    acc_p0_c0 += sdot_vec4xi8(i_oh0_ic2, f_oc0_ic2);
                    acc_p0_c0 += sdot_vec4xi8(i_oh0_ic3, f_oc0_ic3);

                    acc_p1_c0 += sdot_vec4xi8(i_oh1_ic0, f_oc0_ic0);
                    acc_p1_c0 += sdot_vec4xi8(i_oh1_ic1, f_oc0_ic1);
                    acc_p1_c0 += sdot_vec4xi8(i_oh1_ic2, f_oc0_ic2);
                    acc_p1_c0 += sdot_vec4xi8(i_oh1_ic3, f_oc0_ic3);

                    acc_p2_c0 += sdot_vec4xi8(i_oh2_ic0, f_oc0_ic0);
                    acc_p2_c0 += sdot_vec4xi8(i_oh2_ic1, f_oc0_ic1);
                    acc_p2_c0 += sdot_vec4xi8(i_oh2_ic2, f_oc0_ic2);
                    acc_p2_c0 += sdot_vec4xi8(i_oh2_ic3, f_oc0_ic3);

                    acc_p0_c1 += sdot_vec4xi8(i_oh0_ic0, f_oc1_ic0);
                    acc_p0_c1 += sdot_vec4xi8(i_oh0_ic1, f_oc1_ic1);
                    acc_p0_c1 += sdot_vec4xi8(i_oh0_ic2, f_oc1_ic2);
                    acc_p0_c1 += sdot_vec4xi8(i_oh0_ic3, f_oc1_ic3);

                    acc_p1_c1 += sdot_vec4xi8(i_oh1_ic0, f_oc1_ic0);
                    acc_p1_c1 += sdot_vec4xi8(i_oh1_ic1, f_oc1_ic1);
                    acc_p1_c1 += sdot_vec4xi8(i_oh1_ic2, f_oc1_ic2);
                    acc_p1_c1 += sdot_vec4xi8(i_oh1_ic3, f_oc1_ic3);

                    acc_p2_c1 += sdot_vec4xi8(i_oh2_ic0, f_oc1_ic0);
                    acc_p2_c1 += sdot_vec4xi8(i_oh2_ic1, f_oc1_ic1);
                    acc_p2_c1 += sdot_vec4xi8(i_oh2_ic2, f_oc1_ic2);
                    acc_p2_c1 += sdot_vec4xi8(i_oh2_ic3, f_oc1_ic3);

                    in_oh_ic += ICR * ICV;
                    filter_oc_ic += ICR * ICV;
                }  // end ic loop

                out_oc_oh[0 * OC + 0] = acc_p0_c0;
                out_oc_oh[1 * OC + 0] = acc_p1_c0;
                out_oc_oh[2 * OC + 0] = acc_p2_c0;

                out_oc_oh[0 * OC + 1] = acc_p0_c1;
                out_oc_oh[1 * OC + 1] = acc_p1_c1;
                out_oc_oh[2 * OC + 1] = acc_p2_c1;

                out_oc_oh += OHWR * OC;
                in_oh += OHWR * IC;
            }  // end oh loop
            filter_oc += OCR * IC;
            bias_oc += OCR;
            out_oc += OCR;
        }  // end oc loop
    }  // end __effcc_ignore_memory_order
}

void conv2d_layer_5(
    const int8_t *in_base, idx_t in_offset, idx_t in_dim_n, idx_t in_dim_h,
    idx_t in_dim_w, idx_t in_dim_c, idx_t in_sn, idx_t in_sh, idx_t in_sw,
    idx_t in_sc, const int8_t *filter_base, idx_t filter_offset, idx_t fd0,
    idx_t fd1, idx_t fd2, idx_t fd3, idx_t fs0, idx_t fs1, idx_t fs2, idx_t fs3,
    const int32_t *bias_base, idx_t bias_offset, idx_t bias_dim,
    idx_t bias_stride, const int8_t *input_zp_ptr, idx_t input_zp_offset,
    idx_t izp_dim, idx_t izp_stride, const int8_t *weight_zp_ptr,
    idx_t weight_zp_offset, idx_t wzp_dim, idx_t wzp_stride, int32_t *out_base,
    idx_t out_offset, idx_t out_dim_n, idx_t out_dim_h, idx_t out_dim_w,
    idx_t out_dim_c, idx_t out_sn, idx_t out_sh, idx_t out_sw, idx_t out_sc) {
    const int OC = 128;
    const int K = 128;
    int32_t bias_w_zp[OC];
    // assume offsets are zero.
    conv1d_bias_zp(filter_base, bias_base, K, OC, *input_zp_ptr, bias_w_zp);
    opt_conv2d_layer_5(in_base, filter_base, bias_w_zp, out_base);
}

// ============================================================================
// Layer 6: depthwise_conv2d
//   input:  [N][15][1][128] i8
//   weight: [15][1][128][1] i8
//   output: [N][1][1][128]  i32
// ============================================================================
__effcc_rip_exact void opt_dw_conv2d_layer_6(
    const int8_t *in,      // [IH:15, OC:128]:int8
    const int8_t *filter,  // [FH:15, OC:128]:int8
    const int32_t *bias,   // [       OC:128]:int32
    int32_t *out           // [OH: 1, OC:128]:int32
) {
    const int OC = 128;
    const int KH = 15;

    // tile size: vectorize OCV channels at a time
    const int OCV = 4;

    // OH=1, so no OH loop — one call per OC tile
    __effcc_ignore_memory_order {
        const int8_t *in_oc = in;
        const int8_t *filter_oc = filter;
        const int32_t *bias_oc = bias;
        int32_t *out_oc = out;
        for (int oc = 0; oc < OC / OCV; oc++) {
            int32_t b0 = bias_oc[0], b1 = bias_oc[1], b2 = bias_oc[2],
                    b3 = bias_oc[3];

            int32_t o_oc0 = b0, o_oc1 = b1, o_oc2 = b2, o_oc3 = b3;
            const int8_t *in_oc_kh = in_oc;
            const int8_t *filter_oc_kh = filter_oc;
            for (int kh = 0; kh < KH; kh++) {
                word4b_t in_kh0;
                word4b_t f_kh0;

                int32_t mul_oc0, mul_oc1, mul_oc2, mul_oc3;

                in_kh0.u32 = *((const uint32_t *)(in_oc_kh + 0 * OC));

                f_kh0.u32 = *((const uint32_t *)(filter_oc_kh + 0 * OC));

                mul_vec4xi8_retptr_4xi32(in_kh0.u32, f_kh0.u32, &mul_oc0,
                                         &mul_oc1, &mul_oc2, &mul_oc3);
                o_oc0 += mul_oc0;
                o_oc1 += mul_oc1;
                o_oc2 += mul_oc2;
                o_oc3 += mul_oc3;

                in_oc_kh += OC;
                filter_oc_kh += OC;
            }
            out_oc[0] = o_oc0;
            out_oc[1] = o_oc1;
            out_oc[2] = o_oc2;
            out_oc[3] = o_oc3;

            in_oc += OCV;
            filter_oc += OCV;
            bias_oc += OCV;
            out_oc += OCV;
        }
    }
}

void dw_conv2d_layer_6(
    const int8_t *in_base, idx_t in_offset, idx_t in_dim_n, idx_t in_dim_h,
    idx_t in_dim_w, idx_t in_dim_c, idx_t in_sn, idx_t in_sh, idx_t in_sw,
    idx_t in_sc, const int8_t *filter_base, idx_t filter_offset, idx_t fd0,
    idx_t fd1, idx_t fd2, idx_t fd3, idx_t fs0, idx_t fs1, idx_t fs2, idx_t fs3,
    const int32_t *bias_base, idx_t bias_offset, idx_t bias_dim,
    idx_t bias_stride, const int8_t *input_zp_ptr, idx_t input_zp_offset,
    idx_t izp_dim, idx_t izp_stride, const int8_t *weight_zp_ptr,
    idx_t weight_zp_offset, idx_t wzp_dim, idx_t wzp_stride, int32_t *out_base,
    idx_t out_offset, idx_t out_dim_n, idx_t out_dim_h, idx_t out_dim_w,
    idx_t out_dim_c, idx_t out_sn, idx_t out_sh, idx_t out_sw, idx_t out_sc) {
    const int OC = 128;
    const int K = 15;
    int32_t bias_w_zp[OC];
    // assume offsets are zero.
    depthwise_conv1d_bias_zp(filter_base, bias_base, K, OC, *input_zp_ptr,
                             bias_w_zp);
    opt_dw_conv2d_layer_6(in_base, filter_base, bias_w_zp, out_base);
}

// ============================================================================
// Layer 7: conv2d (pointwise 1x1)
//   input:  [N][1][1][128] i8
//   weight: [32][1][1][128] i8
//   output: [N][1][1][32]  i32
// ============================================================================
__effcc_rip_exact void opt_conv2d_layer_7(
    const int8_t *in,      // [IC:128]:int8
    const int8_t *filter,  // [OC:32, IC:128]:int8
    const int32_t *bias,   // [OC:32]:int32
    int32_t *out           // [OC:32]:int32
) {
    const int IC = 128;
    const int OC = 32;
    // vectorize IC by ICV=4 using sdot, unroll IC by ICR=8
    // IC=128 evenly divisible by ICV*ICR=32
    const int ICV = 4;
    const int ICR = 8;
    const int OCR = 1;

    __effcc_ignore_memory_order {
        const int32_t *bias_oc = bias;
        const int8_t *filter_oc = filter;
        int32_t *out_oc = out;

        for (int oc = 0; oc < OC; oc++) {
            int32_t acc0 = bias_oc[0];
            const int8_t *filter_oc_ic = filter_oc;
            const int8_t *in_ic = in;
            for (int ic = 0; ic < IC / (ICV * ICR); ic++) {
                word4b_t in_v0, in_v1, in_v2, in_v3, in_v4, in_v5, in_v6, in_v7;
                word4b_t filt_v0, filt_v1, filt_v2, filt_v3, filt_v4, filt_v5,
                    filt_v6, filt_v7;

                in_v0.u32 = *((const uint32_t *)(in_ic + 0 * ICV));
                in_v1.u32 = *((const uint32_t *)(in_ic + 1 * ICV));
                in_v2.u32 = *((const uint32_t *)(in_ic + 2 * ICV));
                in_v3.u32 = *((const uint32_t *)(in_ic + 3 * ICV));
                in_v4.u32 = *((const uint32_t *)(in_ic + 4 * ICV));
                in_v5.u32 = *((const uint32_t *)(in_ic + 5 * ICV));
                in_v6.u32 = *((const uint32_t *)(in_ic + 6 * ICV));
                in_v7.u32 = *((const uint32_t *)(in_ic + 7 * ICV));

                filt_v0.u32 = *((const uint32_t *)(filter_oc_ic + 0 * ICV));
                filt_v1.u32 = *((const uint32_t *)(filter_oc_ic + 1 * ICV));
                filt_v2.u32 = *((const uint32_t *)(filter_oc_ic + 2 * ICV));
                filt_v3.u32 = *((const uint32_t *)(filter_oc_ic + 3 * ICV));
                filt_v4.u32 = *((const uint32_t *)(filter_oc_ic + 4 * ICV));
                filt_v5.u32 = *((const uint32_t *)(filter_oc_ic + 5 * ICV));
                filt_v6.u32 = *((const uint32_t *)(filter_oc_ic + 6 * ICV));
                filt_v7.u32 = *((const uint32_t *)(filter_oc_ic + 7 * ICV));

                int32_t s0 = sdot_vec4xi8(in_v0.u32, filt_v0.u32) +
                             sdot_vec4xi8(in_v1.u32, filt_v1.u32);
                int32_t s1 = sdot_vec4xi8(in_v2.u32, filt_v2.u32) +
                             sdot_vec4xi8(in_v3.u32, filt_v3.u32);
                int32_t s2 = sdot_vec4xi8(in_v4.u32, filt_v4.u32) +
                             sdot_vec4xi8(in_v5.u32, filt_v5.u32);
                int32_t s3 = sdot_vec4xi8(in_v6.u32, filt_v6.u32) +
                             sdot_vec4xi8(in_v7.u32, filt_v7.u32);
                acc0 += (s0 + s1) + (s2 + s3);

                in_ic += ICR * ICV;
                filter_oc_ic += ICR * ICV;
            }  // end ic loop
            out_oc[0] = acc0;

            bias_oc += OCR;
            filter_oc += OCR * IC;
            out_oc += OCR;
        }  // end oc tile loop
    }  // end __effcc_ignore_memory_order
}

void conv2d_layer_7(
    const int8_t *in_base, idx_t in_offset, idx_t in_dim_n, idx_t in_dim_h,
    idx_t in_dim_w, idx_t in_dim_c, idx_t in_sn, idx_t in_sh, idx_t in_sw,
    idx_t in_sc, const int8_t *filter_base, idx_t filter_offset, idx_t fd0,
    idx_t fd1, idx_t fd2, idx_t fd3, idx_t fs0, idx_t fs1, idx_t fs2, idx_t fs3,
    const int32_t *bias_base, idx_t bias_offset, idx_t bias_dim,
    idx_t bias_stride, const int8_t *input_zp_ptr, idx_t input_zp_offset,
    idx_t izp_dim, idx_t izp_stride, const int8_t *weight_zp_ptr,
    idx_t weight_zp_offset, idx_t wzp_dim, idx_t wzp_stride, int32_t *out_base,
    idx_t out_offset, idx_t out_dim_n, idx_t out_dim_h, idx_t out_dim_w,
    idx_t out_dim_c, idx_t out_sn, idx_t out_sh, idx_t out_sw, idx_t out_sc) {
    const int OC = 32;
    const int K = 128;
    int32_t bias_w_zp[OC];
    conv1d_bias_zp(filter_base, bias_base, K, OC, *input_zp_ptr, bias_w_zp);
    opt_conv2d_layer_7(in_base, filter_base, bias_w_zp, out_base);
}

// ============================================================================
// Layer 8: conv2d (pointwise 1x1, classification head)
//   input:  [N][1][1][32] i8
//   weight: [3][1][1][32] i8
//   output: [N][1][1][3]  i32
// ============================================================================
__effcc_rip_exact void opt_conv2d_layer_8(
    const int8_t *in,      // [IC:32]:int8
    const int8_t *filter,  // [OC:3, IC:32]:int8
    const int32_t *bias,   // [OC:3]:int32
    int32_t *out           // [OC:3]:int32
) {
    // OC=3 fully unrolled; IC=32 fully unrolled via ICV=4, ICR=8 (1 iteration)
    const int IC = 32;
    const int ICV = 4;

    __effcc_ignore_memory_order {
        int32_t acc0 = bias[0], acc1 = bias[1], acc2 = bias[2];

        word4b_t in_v0, in_v1, in_v2, in_v3, in_v4, in_v5, in_v6, in_v7;
        word4b_t filt0_v0, filt0_v1, filt0_v2, filt0_v3, filt0_v4, filt0_v5,
            filt0_v6, filt0_v7;
        word4b_t filt1_v0, filt1_v1, filt1_v2, filt1_v3, filt1_v4, filt1_v5,
            filt1_v6, filt1_v7;
        word4b_t filt2_v0, filt2_v1, filt2_v2, filt2_v3, filt2_v4, filt2_v5,
            filt2_v6, filt2_v7;

        // TODO: this would ahve to be just IC
        in_v0.u32 = *((const uint32_t *)(in + 0 * ICV));
        in_v1.u32 = *((const uint32_t *)(in + 1 * ICV));
        in_v2.u32 = *((const uint32_t *)(in + 2 * ICV));
        in_v3.u32 = *((const uint32_t *)(in + 3 * ICV));
        in_v4.u32 = *((const uint32_t *)(in + 4 * ICV));
        in_v5.u32 = *((const uint32_t *)(in + 5 * ICV));
        in_v6.u32 = *((const uint32_t *)(in + 6 * ICV));
        in_v7.u32 = *((const uint32_t *)(in + 7 * ICV));

        filt0_v0.u32 = *((const uint32_t *)(filter + 0 * IC + 0 * ICV));
        filt0_v1.u32 = *((const uint32_t *)(filter + 0 * IC + 1 * ICV));
        filt0_v2.u32 = *((const uint32_t *)(filter + 0 * IC + 2 * ICV));
        filt0_v3.u32 = *((const uint32_t *)(filter + 0 * IC + 3 * ICV));
        filt0_v4.u32 = *((const uint32_t *)(filter + 0 * IC + 4 * ICV));
        filt0_v5.u32 = *((const uint32_t *)(filter + 0 * IC + 5 * ICV));
        filt0_v6.u32 = *((const uint32_t *)(filter + 0 * IC + 6 * ICV));
        filt0_v7.u32 = *((const uint32_t *)(filter + 0 * IC + 7 * ICV));

        filt1_v0.u32 = *((const uint32_t *)(filter + 1 * IC + 0 * ICV));
        filt1_v1.u32 = *((const uint32_t *)(filter + 1 * IC + 1 * ICV));
        filt1_v2.u32 = *((const uint32_t *)(filter + 1 * IC + 2 * ICV));
        filt1_v3.u32 = *((const uint32_t *)(filter + 1 * IC + 3 * ICV));
        filt1_v4.u32 = *((const uint32_t *)(filter + 1 * IC + 4 * ICV));
        filt1_v5.u32 = *((const uint32_t *)(filter + 1 * IC + 5 * ICV));
        filt1_v6.u32 = *((const uint32_t *)(filter + 1 * IC + 6 * ICV));
        filt1_v7.u32 = *((const uint32_t *)(filter + 1 * IC + 7 * ICV));

        filt2_v0.u32 = *((const uint32_t *)(filter + 2 * IC + 0 * ICV));
        filt2_v1.u32 = *((const uint32_t *)(filter + 2 * IC + 1 * ICV));
        filt2_v2.u32 = *((const uint32_t *)(filter + 2 * IC + 2 * ICV));
        filt2_v3.u32 = *((const uint32_t *)(filter + 2 * IC + 3 * ICV));
        filt2_v4.u32 = *((const uint32_t *)(filter + 2 * IC + 4 * ICV));
        filt2_v5.u32 = *((const uint32_t *)(filter + 2 * IC + 5 * ICV));
        filt2_v6.u32 = *((const uint32_t *)(filter + 2 * IC + 6 * ICV));
        filt2_v7.u32 = *((const uint32_t *)(filter + 2 * IC + 7 * ICV));

        acc0 += sdot_vec4xi8(in_v0.u32, filt0_v0.u32) +
                sdot_vec4xi8(in_v1.u32, filt0_v1.u32) +
                sdot_vec4xi8(in_v2.u32, filt0_v2.u32) +
                sdot_vec4xi8(in_v3.u32, filt0_v3.u32) +
                sdot_vec4xi8(in_v4.u32, filt0_v4.u32) +
                sdot_vec4xi8(in_v5.u32, filt0_v5.u32) +
                sdot_vec4xi8(in_v6.u32, filt0_v6.u32) +
                sdot_vec4xi8(in_v7.u32, filt0_v7.u32);
        acc1 += sdot_vec4xi8(in_v0.u32, filt1_v0.u32) +
                sdot_vec4xi8(in_v1.u32, filt1_v1.u32) +
                sdot_vec4xi8(in_v2.u32, filt1_v2.u32) +
                sdot_vec4xi8(in_v3.u32, filt1_v3.u32) +
                sdot_vec4xi8(in_v4.u32, filt1_v4.u32) +
                sdot_vec4xi8(in_v5.u32, filt1_v5.u32) +
                sdot_vec4xi8(in_v6.u32, filt1_v6.u32) +
                sdot_vec4xi8(in_v7.u32, filt1_v7.u32);
        acc2 += sdot_vec4xi8(in_v0.u32, filt2_v0.u32) +
                sdot_vec4xi8(in_v1.u32, filt2_v1.u32) +
                sdot_vec4xi8(in_v2.u32, filt2_v2.u32) +
                sdot_vec4xi8(in_v3.u32, filt2_v3.u32) +
                sdot_vec4xi8(in_v4.u32, filt2_v4.u32) +
                sdot_vec4xi8(in_v5.u32, filt2_v5.u32) +
                sdot_vec4xi8(in_v6.u32, filt2_v6.u32) +
                sdot_vec4xi8(in_v7.u32, filt2_v7.u32);

        out[0] = acc0;
        out[1] = acc1;
        out[2] = acc2;
    }  // end __effcc_ignore_memory_order
}

void conv2d_layer_8(
    const int8_t *in_base, idx_t in_offset, idx_t in_dim_n, idx_t in_dim_h,
    idx_t in_dim_w, idx_t in_dim_c, idx_t in_sn, idx_t in_sh, idx_t in_sw,
    idx_t in_sc, const int8_t *filter_base, idx_t filter_offset, idx_t fd0,
    idx_t fd1, idx_t fd2, idx_t fd3, idx_t fs0, idx_t fs1, idx_t fs2, idx_t fs3,
    const int32_t *bias_base, idx_t bias_offset, idx_t bias_dim,
    idx_t bias_stride, const int8_t *input_zp_ptr, idx_t input_zp_offset,
    idx_t izp_dim, idx_t izp_stride, const int8_t *weight_zp_ptr,
    idx_t weight_zp_offset, idx_t wzp_dim, idx_t wzp_stride, int32_t *out_base,
    idx_t out_offset, idx_t out_dim_n, idx_t out_dim_h, idx_t out_dim_w,
    idx_t out_dim_c, idx_t out_sn, idx_t out_sh, idx_t out_sw, idx_t out_sc) {
    const int OC = 3;
    const int K = 32;
    int32_t bias_w_zp[OC];
    // assume offsets are zero.
    conv1d_bias_zp(filter_base, bias_base, K, OC, *input_zp_ptr, bias_w_zp);
    opt_conv2d_layer_8(in_base, filter_base, bias_w_zp, out_base);
}
