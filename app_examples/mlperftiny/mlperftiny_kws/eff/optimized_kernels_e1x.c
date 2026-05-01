// <<AUTOBENCH>> efficient_e1x efficient_e1
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include "effintrinsic.h"

#define COMPUTE_ZERO_POINT_ON_DEVICE 1

typedef uint32_t idx_t;

// Types
typedef int8_t vec2xi8_t __attribute__((ext_vector_type(2)));
typedef int16_t vec2xi16_t __attribute__((ext_vector_type(2)));

typedef int8_t vec4xi8_t __attribute__((ext_vector_type(4)));
typedef int16_t vec4xi16_t __attribute__((ext_vector_type(4)));
typedef int32_t vec4xi32_t __attribute__((ext_vector_type(4)));

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

extern int32_t sdot_vec4xi8(vec4xi8_t a, vec4xi8_t b);
extern uint16_t load_guard_i16(bool valid, uint16_t a, uint16_t b);

extern void mul_vec4xi8_retptr_4xi32(uint32_t a, uint32_t b, int32_t *c0,
                                     int32_t *c1, int32_t *c2, int32_t *c3);

typedef union {
    vec4xi16_t vi16;
    int16_t ai16[4];
    uint64_t u64;
} word8b_t;

typedef union {
    vec4xi8_t vi8;
    vec2xi16_t vi16;
    int8_t ai8[4];
    uint32_t u32;
} word4b_t;

typedef union {
    vec2xi8_t vi8;
    int8_t ai8[2];
    uint16_t u16;
} word2b_t;
// TODO: vectorize these.
// Filter: [OC, FH, FW, IC]
// Bias  : [OC]
__efficient__ void conv2d_bias_zp(const int8_t *filter, idx_t filter_stride_oc,
                                  const int32_t *bias, idx_t bias_dim_oc,
                                  int32_t zp, int32_t *bias_w_zp) {
    __effcc_ignore_memory_order {
        for (idx_t i = 0; i < bias_dim_oc; ++i) {
            int32_t bias_zp = 0;
            __effcc_parallel(16) for (idx_t j = 0; j < filter_stride_oc; j++) {
                bias_zp += (int32_t)filter[j];
            }
            filter += filter_stride_oc;
            bias_w_zp[i] = bias[i] + zp * bias_zp;
        }
    }
}

// Filter: [FH, FW, OC] depthwise (OC==IC)
// Bias  : [OC]
__efficient__ void depthwise_conv2d_bias_zp(const int8_t *filter,
                                            const int32_t *bias, idx_t OC,
                                            int32_t zp, int32_t *bias_w_zp) {
    __effcc_ignore_memory_order {
        const idx_t KH = 3;
        const idx_t KW = 3;
        __effcc_parallel(16) for (idx_t oc = 0; oc < OC; oc++) {
            int32_t weight_sum = 0;
            for (idx_t khw = 0; khw < KH * KW; khw++) {
                weight_sum += (int32_t)filter[khw * OC + oc];
            }
            bias_w_zp[oc] = bias[oc] + zp * weight_sum;
        }
    }
}

// populate a 4 byte word via 1+2+1 bytes
#define PACK_4_FROM_121(a, bc, d)                           \
    ((uint32_t)(uint8_t)a | ((uint32_t)(uint16_t)bc << 8) | \
     ((uint32_t)(uint8_t)d << 24));

// NOTE: Assumes zero-point has been folded into bias.
__effcc_rip_exact void conv_kws_layer_0(const int8_t *restrict in,
                                        const int8_t *restrict filter,
                                        const int32_t *restrict bias,
                                        int32_t *restrict out) {
    // --- Parameters ---
    // input size
    const uint32_t IH = 49;
    const uint32_t IW = 10;
    const uint32_t IC = 1;
    // input zero-point
    const int8_t IZP = 83;

    // filter size
    const uint32_t FH = 10;
    const uint32_t FW = 4;
    const uint32_t OC = 64;

    // stride
    const uint32_t ST = 2;

    // padding
    const uint32_t PT = 4;
    const uint32_t PB = 5;
    const uint32_t PL = 1;
    const uint32_t PR = 1;

    // output size
    const uint32_t OH = (IH + PT + PB - FH) / ST + 1;  // 25
    const uint32_t OW = (IW + PL + PR - FW) / ST + 1;  // 5

    // Tile sizes
    const int FHR = 2;
    const int OCR = 1;

    /*
        Pointer convention:
        + use pointer bumping for each loop
            - changes index * stride arithmetic to += stride
            - assume any constant offsets to pointers is okay
        + append loop iterator name to pointer name
            - e.g. out_oc    is equivalent to out    + oc * <oc stride>
            - e.g. out_oc_oh is equivalent to out_oc + oh * <oh stride>
    */
    __effcc_ignore_memory_order {
        const int8_t *filter_oc = filter;
        const int32_t *bias_oc = bias;
        int32_t *out_oc = out;
        for (int oc = 0; oc < OC; oc += OCR) {
            int in_row = 0;
            const int8_t *in_oh = in;
            int32_t *out_oc_oh = out_oc;
            for (int oh = 0; oh < OH; oh++) {
                // fully unroll OW, FW
                int32_t acc_ow0 = bias_oc[0];
                int32_t acc_ow1 = acc_ow0;
                int32_t acc_ow2 = acc_ow0;
                int32_t acc_ow3 = acc_ow0;
                int32_t acc_ow4 = acc_ow0;

                // FH range that map to valid input; otherwise map to pad
                int fh_start = MAX(((0 + PT) - in_row), 0);
                int fh_end = MIN(((IH + PT) - in_row), FH);
                const int8_t *in_oh_fh = in_oh - PT * IW;
                const int8_t *filter_oc_fh = filter_oc;
                for (int fh = 0; fh < FH; fh += FHR) {
                    /*
                        Variable naming convention:
                            in_fh<FHR unroll index>_<FW vectorized indices>
                                - Inputs to 4-wide vector dot product.
                            filt_fh<FHR unroll index>
                                - Weights to 4-wide vector dot product.

                        Vectorizes FW=4 dimension of filter, which gets mapped
                       to 4 element windows from the input width. (e.g. elements
                       1-4 in the width would be _1234).

                        Use z to correspond to padded input (initialized to
                       IZP).

                        Different output elements are at a stride of 2,
                        (e.g. skip _2345, use _3456).

                        We compute FHR=2 different filter rows in each iteration
                        of the FH loop (e.g. in_fh0_* and in_fh1_*),
                        which get accumulated to the same output element in the
                        output width (e.g. acc_ow<OW index>).
                        */
                    word4b_t in_fh0_z012, in_fh0_1234, in_fh0_3456, in_fh0_5678,
                        in_fh0_789z, in_fh1_z012, in_fh1_1234, in_fh1_3456,
                        in_fh1_5678, in_fh1_789z, f_oc0_fh0, f_oc0_fh1;

                    // No packing for filter required since data is 4B aligned
                    // and memory layout matches instruction layout.
                    f_oc0_fh0.u32 =
                        *((uint32_t *)(filter_oc_fh + 0 * FH * FW + 0 * FW));
                    f_oc0_fh1.u32 =
                        *((uint32_t *)(filter_oc_fh + 0 * FH * FW + 1 * FW));

                    bool row0_valid = (fh_start <= fh + 0 && fh + 0 < fh_end);
                    bool row1_valid = (fh_start <= fh + 1 && fh + 1 < fh_end);

                    word2b_t in_fh0_01, in_fh0_23, in_fh0_45, in_fh0_67,
                        in_fh0_89, in_fh1_01, in_fh1_23, in_fh1_45, in_fh1_67,
                        in_fh1_89, zp_v;

                    // Load 2 elements at a time via i16 to reduce bank
                    // conflicts when accessing different i8s from the same
                    // word.
                    // Vectorizing by 4 seemed to cause alignment issues
                    // since (IW = 10) % 4 != 0
                    zp_v.vi8[0] = IZP;
                    zp_v.vi8[1] = IZP;
                    in_fh0_01.u16 = load_guard_i16(
                        row0_valid, *((uint16_t *)(in_oh_fh + 0 * IW + 0)),
                        zp_v.u16);
                    in_fh0_23.u16 = load_guard_i16(
                        row0_valid, *((uint16_t *)(in_oh_fh + 0 * IW + 2)),
                        zp_v.u16);
                    in_fh0_45.u16 = load_guard_i16(
                        row0_valid, *((uint16_t *)(in_oh_fh + 0 * IW + 4)),
                        zp_v.u16);
                    in_fh0_67.u16 = load_guard_i16(
                        row0_valid, *((uint16_t *)(in_oh_fh + 0 * IW + 6)),
                        zp_v.u16);
                    in_fh0_89.u16 = load_guard_i16(
                        row0_valid, *((uint16_t *)(in_oh_fh + 0 * IW + 8)),
                        zp_v.u16);

                    in_fh1_01.u16 = load_guard_i16(
                        row1_valid, *((uint16_t *)(in_oh_fh + 1 * IW + 0)),
                        zp_v.u16);
                    in_fh1_23.u16 = load_guard_i16(
                        row1_valid, *((uint16_t *)(in_oh_fh + 1 * IW + 2)),
                        zp_v.u16);
                    in_fh1_45.u16 = load_guard_i16(
                        row1_valid, *((uint16_t *)(in_oh_fh + 1 * IW + 4)),
                        zp_v.u16);
                    in_fh1_67.u16 = load_guard_i16(
                        row1_valid, *((uint16_t *)(in_oh_fh + 1 * IW + 6)),
                        zp_v.u16);
                    in_fh1_89.u16 = load_guard_i16(
                        row1_valid, *((uint16_t *)(in_oh_fh + 1 * IW + 8)),
                        zp_v.u16);

                    // Repack input pairs of elements into the 4-byte words
                    // used for vector dot product. Allows us to reuse loads
                    // rather than loading the same elements redundantly.
                    in_fh0_z012.u32 =
                        PACK_4_FROM_121(IZP, in_fh0_01.u16, in_fh0_23.ai8[0]);
                    in_fh0_1234.u32 = PACK_4_FROM_121(
                        in_fh0_01.ai8[1], in_fh0_23.u16, in_fh0_45.ai8[0]);
                    in_fh0_3456.u32 = PACK_4_FROM_121(
                        in_fh0_23.ai8[1], in_fh0_45.u16, in_fh0_67.ai8[0]);
                    in_fh0_5678.u32 = PACK_4_FROM_121(
                        in_fh0_45.ai8[1], in_fh0_67.u16, in_fh0_89.ai8[0]);
                    in_fh0_789z.u32 =
                        PACK_4_FROM_121(in_fh0_67.ai8[1], in_fh0_89.u16, IZP);

                    in_fh1_z012.u32 =
                        PACK_4_FROM_121(IZP, in_fh1_01.u16, in_fh1_23.ai8[0]);
                    in_fh1_1234.u32 = PACK_4_FROM_121(
                        in_fh1_01.ai8[1], in_fh1_23.u16, in_fh1_45.ai8[0]);
                    in_fh1_3456.u32 = PACK_4_FROM_121(
                        in_fh1_23.ai8[1], in_fh1_45.u16, in_fh1_67.ai8[0]);
                    in_fh1_5678.u32 = PACK_4_FROM_121(
                        in_fh1_45.ai8[1], in_fh1_67.u16, in_fh1_89.ai8[0]);
                    in_fh1_789z.u32 =
                        PACK_4_FROM_121(in_fh1_67.ai8[1], in_fh1_89.u16, IZP);

                    acc_ow0 += sdot_vec4xi8(in_fh0_z012.vi8, f_oc0_fh0.vi8);
                    acc_ow1 += sdot_vec4xi8(in_fh0_1234.vi8, f_oc0_fh0.vi8);
                    acc_ow2 += sdot_vec4xi8(in_fh0_3456.vi8, f_oc0_fh0.vi8);
                    acc_ow3 += sdot_vec4xi8(in_fh0_5678.vi8, f_oc0_fh0.vi8);
                    acc_ow4 += sdot_vec4xi8(in_fh0_789z.vi8, f_oc0_fh0.vi8);

                    acc_ow0 += sdot_vec4xi8(in_fh1_z012.vi8, f_oc0_fh1.vi8);
                    acc_ow1 += sdot_vec4xi8(in_fh1_1234.vi8, f_oc0_fh1.vi8);
                    acc_ow2 += sdot_vec4xi8(in_fh1_3456.vi8, f_oc0_fh1.vi8);
                    acc_ow3 += sdot_vec4xi8(in_fh1_5678.vi8, f_oc0_fh1.vi8);
                    acc_ow4 += sdot_vec4xi8(in_fh1_789z.vi8, f_oc0_fh1.vi8);

                    in_oh_fh += FHR * IW;
                    filter_oc_fh += FHR * FW;
                }
                out_oc_oh[0 * OC + 0] = acc_ow0;
                out_oc_oh[1 * OC + 0] = acc_ow1;
                out_oc_oh[2 * OC + 0] = acc_ow2;
                out_oc_oh[3 * OC + 0] = acc_ow3;
                out_oc_oh[4 * OC + 0] = acc_ow4;

                in_row += ST;
                in_oh += ST * IW;
                out_oc_oh += OW * OC;
            }  // end oh loop
            filter_oc += OCR * FH * FW;
            bias_oc += OCR;
            out_oc += OCR;
        }  // end oc loop
    }  // end __effcc_ignore_memory_order
}

void conv2d_kws_layer_0(
    // ==== Input (i8, NHWC) ====
    const int8_t *in_base, idx_t in_offset, idx_t in_dim_n, idx_t in_dim_h,
    idx_t in_dim_w, idx_t in_dim_c, idx_t in_stride_n, idx_t in_stride_h,
    idx_t in_stride_w, idx_t in_stride_c,

    // ==== Filter (i8, [OC,KH,KW,IC]) ====
    const int8_t *filter_base, idx_t filter_offset, idx_t filter_dim_oc,
    idx_t filter_dim_kh, idx_t filter_dim_kw, idx_t filter_dim_ic,
    idx_t filter_stride_oc, idx_t filter_stride_h, idx_t filter_stride_w,
    idx_t filter_stride_ic,

    // ==== Bias (i32, length = OC) ====
    const int32_t *bias_base, idx_t bias_offset, idx_t bias_dim_oc,
    idx_t bias_stride_1,

    // ==== input_zp (i8, length = 1) ====
    const int8_t *input_zp_ptr, idx_t input_zp_offset, idx_t input_zp_dim_1,
    idx_t input_zp_stride_1,

    // ==== weight_zp (i8, length = 1) ====
    const int8_t *weight_zp_ptr, idx_t weight_zp_offset, idx_t weight_zp_dim_1,
    idx_t weight_zp_stride_1,

    // ==== Output (i32, NHWC) ====
    int32_t *out_base, idx_t out_offset, idx_t out_dim_n, idx_t out_dim_h,
    idx_t out_dim_w, idx_t out_dim_c, idx_t out_stride_n, idx_t out_stride_h,
    idx_t out_stride_w, idx_t out_stride_c) {
#if COMPUTE_ZERO_POINT_ON_DEVICE
    // zero-point folding into bias
    const int OC = 64;
    int32_t bias_w_zp[OC];
    conv2d_bias_zp(filter_base, filter_stride_oc, bias_base, filter_dim_oc,
                   -input_zp_ptr[0], bias_w_zp);

    conv_kws_layer_0(in_base, filter_base, bias_w_zp, out_base);
#else
    conv_kws_layer_0(in_base, filter_base, bias_base, out_base);
#endif
}

// NOTE: Assumes zero-point has been folded into bias.
__effcc_rip_exact void depthwise_conv_kws_layer_1(
    const int8_t *restrict in,      // [IH][IW][IC] NHWC layout (assuming N==1)
    const int8_t *restrict filter,  // [KH][KW][IC=OC][depth_mult=1]
    const int32_t *restrict bias,   // [OC] (OC == IC)
    int32_t *restrict out           // [OH][OW][OC], NHWC
) {
    // --- Parameters ---
    // input size
    const uint32_t IH = 25;
    const uint32_t IW = 5;
    const uint32_t IC = 64;
    // input zero-point
    const int8_t IZP = -128;

    // filter size
    const uint32_t FH = 3;
    const uint32_t FW = 3;
    const uint32_t OC = 64;

    // stride
    const uint32_t ST = 1;

    // padding
    const uint32_t PT = 1;
    const uint32_t PB = 1;
    const uint32_t PL = 1;
    const uint32_t PR = 1;

    // output size
    const uint32_t OH = (IH + PT + PB - FH) / ST + 1;  // 25
    const uint32_t OW = (IW + PL + PR - FW) / ST + 1;  // 5

    // tile size, vectorize loads of int8_t values from inputs, weights.
    const uint32_t OCR = 4;

    __effcc_ignore_memory_order {
        const int32_t *bias_oc = bias;
        int32_t *out_oc = out;
        const int8_t *in_oc = in;
        const int8_t *filter_oc = filter;
        for (int oc = 0; oc < OC / OCR; oc++) {
            int32_t b_oc0 = bias_oc[0], b_oc1 = bias_oc[1], b_oc2 = bias_oc[2],
                    b_oc3 = bias_oc[3];
            int32_t *out_oc_oh = out_oc;
            const int8_t *in_oc_oh = in_oc;
            for (int oh = 0; oh < OH; oh++) {
                // fully unroll OW, FW
                int32_t o_ow0_oc0 = b_oc0, o_ow0_oc1 = b_oc1, o_ow0_oc2 = b_oc2,
                        o_ow0_oc3 = b_oc3;
                int32_t o_ow1_oc0 = b_oc0, o_ow1_oc1 = b_oc1, o_ow1_oc2 = b_oc2,
                        o_ow1_oc3 = b_oc3;
                int32_t o_ow2_oc0 = b_oc0, o_ow2_oc1 = b_oc1, o_ow2_oc2 = b_oc2,
                        o_ow2_oc3 = b_oc3;
                int32_t o_ow3_oc0 = b_oc0, o_ow3_oc1 = b_oc1, o_ow3_oc2 = b_oc2,
                        o_ow3_oc3 = b_oc3;
                int32_t o_ow4_oc0 = b_oc0, o_ow4_oc1 = b_oc1, o_ow4_oc2 = b_oc2,
                        o_ow4_oc3 = b_oc3;

                int fh_start = (PT - oh) > 0 ? (PT - oh) : 0;
                int fh_end = ((IH + PT) - oh) < FH ? ((IH + PT) - oh) : FH;
                const int8_t *in_oc_oh_fh = in_oc_oh - PT * IW * IC;
                const int8_t *filter_oc_fh = filter_oc;
                for (int fh = 0; fh < FH; fh++) {
                    bool valid = fh_start <= fh && fh < fh_end;
                    word4b_t in_iw0_v, in_iw1_v, in_iw2_v, in_iw3_v, in_iw4_v,
                        filt_fw0, filt_fw1, filt_fw2, zp_v, win0, win1, win2,
                        win3, win4, win5, win6;
                    zp_v.ai8[0] = IZP;
                    zp_v.ai8[1] = IZP;
                    zp_v.ai8[2] = IZP;
                    zp_v.ai8[3] = IZP;

                    in_iw0_v.u32 = valid ? *((uint32_t *)(in_oc_oh_fh + 0 * IC))
                                         : zp_v.u32;
                    in_iw1_v.u32 = valid ? *((uint32_t *)(in_oc_oh_fh + 1 * IC))
                                         : zp_v.u32;
                    in_iw2_v.u32 = valid ? *((uint32_t *)(in_oc_oh_fh + 2 * IC))
                                         : zp_v.u32;
                    in_iw3_v.u32 = valid ? *((uint32_t *)(in_oc_oh_fh + 3 * IC))
                                         : zp_v.u32;
                    in_iw4_v.u32 = valid ? *((uint32_t *)(in_oc_oh_fh + 4 * IC))
                                         : zp_v.u32;

                    // Sliding window: win0=left_pad, win1..win5=iw0..iw4,
                    // win6=right_pad
                    win0 = zp_v;
                    win1 = in_iw0_v;
                    win2 = in_iw1_v;
                    win3 = in_iw2_v;
                    win4 = in_iw3_v;
                    win5 = in_iw4_v;
                    win6 = zp_v;

                    filt_fw0.u32 = *((uint32_t *)(filter_oc_fh + 0 * OC));
                    filt_fw1.u32 = *((uint32_t *)(filter_oc_fh + 1 * OC));
                    filt_fw2.u32 = *((uint32_t *)(filter_oc_fh + 2 * OC));

                    // vectorize OCR, unroll IW
                    for (int fw = 0; fw < FW; fw++) {
                        int32_t ow00, ow01, ow02, ow03, ow10, ow11, ow12, ow13,
                            ow20, ow21, ow22, ow23, ow30, ow31, ow32, ow33,
                            ow40, ow41, ow42, ow43;

                        mul_vec4xi8_retptr_4xi32(win0.u32, filt_fw0.u32, &ow00,
                                                 &ow01, &ow02, &ow03);
                        o_ow0_oc0 += ow00;
                        o_ow0_oc1 += ow01;
                        o_ow0_oc2 += ow02;
                        o_ow0_oc3 += ow03;

                        mul_vec4xi8_retptr_4xi32(win1.u32, filt_fw0.u32, &ow10,
                                                 &ow11, &ow12, &ow13);
                        o_ow1_oc0 += ow10;
                        o_ow1_oc1 += ow11;
                        o_ow1_oc2 += ow12;
                        o_ow1_oc3 += ow13;

                        mul_vec4xi8_retptr_4xi32(win2.u32, filt_fw0.u32, &ow20,
                                                 &ow21, &ow22, &ow23);
                        o_ow2_oc0 += ow20;
                        o_ow2_oc1 += ow21;
                        o_ow2_oc2 += ow22;
                        o_ow2_oc3 += ow23;

                        mul_vec4xi8_retptr_4xi32(win3.u32, filt_fw0.u32, &ow30,
                                                 &ow31, &ow32, &ow33);
                        o_ow3_oc0 += ow30;
                        o_ow3_oc1 += ow31;
                        o_ow3_oc2 += ow32;
                        o_ow3_oc3 += ow33;

                        mul_vec4xi8_retptr_4xi32(win4.u32, filt_fw0.u32, &ow40,
                                                 &ow41, &ow42, &ow43);
                        o_ow4_oc0 += ow40;
                        o_ow4_oc1 += ow41;
                        o_ow4_oc2 += ow42;
                        o_ow4_oc3 += ow43;

                        win0 = win1;
                        win1 = win2;
                        win2 = win3;
                        win3 = win4;
                        win4 = win5;
                        win5 = win6;
                        filt_fw0 = filt_fw1;
                        filt_fw1 = filt_fw2;
                    }
                    filter_oc_fh += FW * OC;
                    in_oc_oh_fh += IW * IC;
                }  // end fh loop
                out_oc_oh[0 * OC + 0] = o_ow0_oc0;
                out_oc_oh[0 * OC + 1] = o_ow0_oc1;
                out_oc_oh[0 * OC + 2] = o_ow0_oc2;
                out_oc_oh[0 * OC + 3] = o_ow0_oc3;

                out_oc_oh[1 * OC + 0] = o_ow1_oc0;
                out_oc_oh[1 * OC + 1] = o_ow1_oc1;
                out_oc_oh[1 * OC + 2] = o_ow1_oc2;
                out_oc_oh[1 * OC + 3] = o_ow1_oc3;

                out_oc_oh[2 * OC + 0] = o_ow2_oc0;
                out_oc_oh[2 * OC + 1] = o_ow2_oc1;
                out_oc_oh[2 * OC + 2] = o_ow2_oc2;
                out_oc_oh[2 * OC + 3] = o_ow2_oc3;

                out_oc_oh[3 * OC + 0] = o_ow3_oc0;
                out_oc_oh[3 * OC + 1] = o_ow3_oc1;
                out_oc_oh[3 * OC + 2] = o_ow3_oc2;
                out_oc_oh[3 * OC + 3] = o_ow3_oc3;

                out_oc_oh[4 * OC + 0] = o_ow4_oc0;
                out_oc_oh[4 * OC + 1] = o_ow4_oc1;
                out_oc_oh[4 * OC + 2] = o_ow4_oc2;
                out_oc_oh[4 * OC + 3] = o_ow4_oc3;

                out_oc_oh += OW * OC;
                in_oc_oh += IW * IC;
            }  // end oh loop
            in_oc += OCR;
            bias_oc += OCR;
            out_oc += OCR;
            filter_oc += OCR;
        }  // end oc loop
    }  // end __effcc_ignore_memory_order
}

void depthwise_conv2d_kws_layer_1(
    // ==== Input (i8, NHWC) ====
    const int8_t *in_base, idx_t in_offset, idx_t in_dim_n, idx_t in_dim_h,
    idx_t in_dim_w, idx_t in_dim_c, idx_t in_stride_n, idx_t in_stride_h,
    idx_t in_stride_w, idx_t in_stride_c,

    // ==== Filter (i8, [OC,KH,KW,IC]) ====
    const int8_t *filter_base, idx_t filter_offset, idx_t filter_dim_oc,
    idx_t filter_dim_kh, idx_t filter_dim_kw, idx_t filter_dim_ic,
    idx_t filter_stride_oc, idx_t filter_stride_h, idx_t filter_stride_w,
    idx_t filter_stride_ic,

    // ==== Bias (i32, length = OC) ====
    const int32_t *bias_base, idx_t bias_offset, idx_t bias_dim_oc,
    idx_t bias_stride_1,

    // ==== input_zp (i8, length = 1) ====
    const int8_t *input_zp_ptr, idx_t input_zp_offset, idx_t input_zp_dim_1,
    idx_t input_zp_stride_1,

    // ==== weight_zp (i8, length = 1) ====
    const int8_t *weight_zp_ptr, idx_t weight_zp_offset, idx_t weight_zp_dim_1,
    idx_t weight_zp_stride_1,

    // ==== Output (i32, NHWC) ====
    int32_t *out_base, idx_t out_offset, idx_t out_dim_n, idx_t out_dim_h,
    idx_t out_dim_w, idx_t out_dim_c, idx_t out_stride_n, idx_t out_stride_h,
    idx_t out_stride_w, idx_t out_stride_c) {
#if COMPUTE_ZERO_POINT_ON_DEVICE
    // zero-point folding into bias
    const int OC = 64;
    int32_t bias_w_zp[OC];
    depthwise_conv2d_bias_zp(filter_base + filter_offset, bias_base, OC,
                             -(int32_t)input_zp_ptr[0], bias_w_zp);
    depthwise_conv_kws_layer_1(in_base, filter_base + filter_offset, bias_w_zp,
                               out_base);
#else
    depthwise_conv_kws_layer_1(in_base, filter_base + filter_offset, bias_base,
                               out_base);
#endif
}

// NOTE: Assumes zero-point has been folded into bias.
__effcc_rip_exact void conv2d_1x1(
    const int8_t *restrict in,      // [IH][IW][IC] NHWC layout (assuming N==1)
    const int8_t *restrict filter,  // [OC][KH][KW][IC]
    const int32_t *restrict bias,   // [OC]
    int32_t *restrict out           // [OH][OW][OC], NHWC
) {
    const uint32_t IH = 25;
    const uint32_t IW = 5;
    const uint32_t IC = 64;

    // No padding so no zero point needed in this kernel
    // const int8_t IZP = -128;

    const uint32_t FH = 1;
    const uint32_t FW = 1;
    const uint32_t OC = 64;

    // merge OH, OW into OHW
    const uint32_t OH = 25;
    const uint32_t OW = 5;
    // logically merge these dimensions
    const uint32_t OHW = OH * OW;

    // vectorize ICV via dot product
    // unroll OHWR, OCR
    const uint32_t ICV = 4;
    const uint32_t OHWR = 3;
    const uint32_t OCR = 4;

    // in/filter/bias are read-only
    // out is write-only (each output is fully computed before writing)
    __effcc_ignore_memory_order {
        int32_t *out_oc = out;
        const int32_t *bias_oc = bias;
        const int8_t *filter_oc = filter;
        for (int oc = 0; oc < (OC / OCR); oc++) {
            const int8_t *in_ohw = in;
            int32_t *out_oc_ohw = out_oc;
            for (int ohw = 0; ohw < ((OHW + OHWR - 1) / OHWR); ohw++) {
                // Accumulators: acc_p{position in ohw tile}_c{channel in oc
                // tile}
                int32_t acc_p0_c0 = bias_oc[0];
                int32_t acc_p1_c0 = acc_p0_c0;
                int32_t acc_p2_c0 = acc_p0_c0;

                int32_t acc_p0_c1 = bias_oc[1];
                int32_t acc_p1_c1 = acc_p0_c1;
                int32_t acc_p2_c1 = acc_p0_c1;

                int32_t acc_p0_c2 = bias_oc[2];
                int32_t acc_p1_c2 = acc_p0_c2;
                int32_t acc_p2_c2 = acc_p0_c2;

                int32_t acc_p0_c3 = bias_oc[3];
                int32_t acc_p1_c3 = acc_p0_c3;
                int32_t acc_p2_c3 = acc_p0_c3;

                const int8_t *filter_oc_ic = filter_oc;
                const int8_t *in_ohw_ic = in_ohw;
                for (int ic = 0; ic < (IC / ICV); ic++) {
                    word4b_t in_p0_v, in_p1_v, in_p2_v, filt_c0_v, filt_c1_v,
                        filt_c2_v, filt_c3_v;
                    in_p0_v.u32 = *((uint32_t *)(in_ohw_ic + 0 * IC));
                    in_p1_v.u32 = *((uint32_t *)(in_ohw_ic + 1 * IC));
                    // know that 125 % 3 is 2 < 3, so only check last row
                    in_p2_v.u32 = __builtin_expect(ohw < ((OHW) / OHWR), 1)
                                      ? *((uint32_t *)(in_ohw_ic + 2 * IC))
                                      : 0;

                    filt_c0_v.u32 = *((uint32_t *)(filter_oc_ic + 0 * IC));
                    filt_c1_v.u32 = *((uint32_t *)(filter_oc_ic + 1 * IC));
                    filt_c2_v.u32 = *((uint32_t *)(filter_oc_ic + 2 * IC));
                    filt_c3_v.u32 = *((uint32_t *)(filter_oc_ic + 3 * IC));

                    acc_p0_c0 += sdot_vec4xi8(in_p0_v.vi8, filt_c0_v.vi8);
                    acc_p0_c1 += sdot_vec4xi8(in_p0_v.vi8, filt_c1_v.vi8);
                    acc_p0_c2 += sdot_vec4xi8(in_p0_v.vi8, filt_c2_v.vi8);
                    acc_p0_c3 += sdot_vec4xi8(in_p0_v.vi8, filt_c3_v.vi8);

                    acc_p1_c0 += sdot_vec4xi8(in_p1_v.vi8, filt_c0_v.vi8);
                    acc_p1_c1 += sdot_vec4xi8(in_p1_v.vi8, filt_c1_v.vi8);
                    acc_p1_c2 += sdot_vec4xi8(in_p1_v.vi8, filt_c2_v.vi8);
                    acc_p1_c3 += sdot_vec4xi8(in_p1_v.vi8, filt_c3_v.vi8);

                    acc_p2_c0 += sdot_vec4xi8(in_p2_v.vi8, filt_c0_v.vi8);
                    acc_p2_c1 += sdot_vec4xi8(in_p2_v.vi8, filt_c1_v.vi8);
                    acc_p2_c2 += sdot_vec4xi8(in_p2_v.vi8, filt_c2_v.vi8);
                    acc_p2_c3 += sdot_vec4xi8(in_p2_v.vi8, filt_c3_v.vi8);

                    in_ohw_ic += ICV;
                    filter_oc_ic += ICV;
                }  // end ic loop
                out_oc_ohw[0 * OC + 0] = acc_p0_c0;
                out_oc_ohw[0 * OC + 1] = acc_p0_c1;
                out_oc_ohw[0 * OC + 2] = acc_p0_c2;
                out_oc_ohw[0 * OC + 3] = acc_p0_c3;

                out_oc_ohw[1 * OC + 0] = acc_p1_c0;
                out_oc_ohw[1 * OC + 1] = acc_p1_c1;
                out_oc_ohw[1 * OC + 2] = acc_p1_c2;
                out_oc_ohw[1 * OC + 3] = acc_p1_c3;

                if (__builtin_expect(ohw < (OHW / OHWR), 1)) {
                    out_oc_ohw[2 * OC + 0] = acc_p2_c0;
                    out_oc_ohw[2 * OC + 1] = acc_p2_c1;
                    out_oc_ohw[2 * OC + 2] = acc_p2_c2;
                    out_oc_ohw[2 * OC + 3] = acc_p2_c3;
                }
                out_oc_ohw += OHWR * OC;
                in_ohw += OHWR * IC;
            }  // end ohw loop
            filter_oc += OCR * IC;
            bias_oc += OCR;
            out_oc += OCR;
        }  // end oc loop
    }
}

void conv2d_kws_layer_2(
    // ==== Input (i8, NHWC) ====
    const int8_t *in_base, idx_t in_offset, idx_t in_dim_n, idx_t in_dim_h,
    idx_t in_dim_w, idx_t in_dim_c, idx_t in_stride_n, idx_t in_stride_h,
    idx_t in_stride_w, idx_t in_stride_c,

    // ==== Filter (i8, [OC,KH,KW,IC]) ====
    const int8_t *filter_base, idx_t filter_offset, idx_t filter_dim_oc,
    idx_t filter_dim_kh, idx_t filter_dim_kw, idx_t filter_dim_ic,
    idx_t filter_stride_oc, idx_t filter_stride_h, idx_t filter_stride_w,
    idx_t filter_stride_ic,

    // ==== Bias (i32, length = OC) ====
    const int32_t *bias_base, idx_t bias_offset, idx_t bias_dim_oc,
    idx_t bias_stride_1,

    // ==== input_zp (i8, length = 1) ====
    const int8_t *input_zp_ptr, idx_t input_zp_offset, idx_t input_zp_dim_1,
    idx_t input_zp_stride_1,

    // ==== weight_zp (i8, length = 1) ====
    const int8_t *weight_zp_ptr, idx_t weight_zp_offset, idx_t weight_zp_dim_1,
    idx_t weight_zp_stride_1,

    // ==== Output (i32, NHWC) ====
    int32_t *out_base, idx_t out_offset, idx_t out_dim_n, idx_t out_dim_h,
    idx_t out_dim_w, idx_t out_dim_c, idx_t out_stride_n, idx_t out_stride_h,
    idx_t out_stride_w, idx_t out_stride_c) {
#if COMPUTE_ZERO_POINT_ON_DEVICE
    // zero-point folding into bias
    const int OC = 64;
    int32_t bias_w_zp[OC];
    conv2d_bias_zp(filter_base, filter_stride_oc, bias_base, filter_dim_oc,
                   -input_zp_ptr[0], bias_w_zp);

    conv2d_1x1(in_base, filter_base, bias_w_zp, out_base);
#else
    conv2d_1x1(in_base, filter_base, bias_base, out_base);
#endif
}
