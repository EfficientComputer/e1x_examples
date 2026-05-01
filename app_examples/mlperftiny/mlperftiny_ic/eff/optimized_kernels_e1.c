// <<AUTOBENCH>> efficient_e1
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include <stdlib.h>

#include "effintrinsic.h"

typedef uint32_t idx_t;

#define likely(x) __builtin_expect(!!(x), 1)

typedef int8_t vec4xi8_t __attribute__((ext_vector_type(4)));

extern int32_t sdot_vec4xi8(vec4xi8_t a, vec4xi8_t b);
extern int32_t add_tree_8(int32_t a, int32_t b, int32_t c, int32_t d, int32_t e,
                          int32_t f, int32_t g, int32_t h);
extern int32_t add_tree_4(int32_t a, int32_t b, int32_t c, int32_t d);
extern vec4xi8_t load_guard(bool cond, vec4xi8_t ld, vec4xi8_t alt);

const idx_t OC_unroll_2 = 2;
const idx_t IC_increment_4 = 4;
const idx_t IC_increment_16 = 16;
const idx_t CH_increment_4 = 4;

/*
 * Same loop order as original (oc → ic → h → w), but 2 OC at once.
 * Minimal diff from _w_const: same input walking, same sliding window,
 * just doubled filter taps and output writes.
 */
__efficient__ void conv2d_stride_1_pad_1_impl_w_const_oc_unroll2(
    const int8_t *input, idx_t IW, const int8_t *filter, idx_t OC, idx_t IC,
    int8_t zp, const int32_t *bias, int32_t *output, idx_t OH, idx_t OW,
    const idx_t IC_div4, const idx_t IW_IC_div4, const idx_t IW_IC_div4_2,
    vec4xi8_t zp_vec, const int8_t *input_base, const idx_t IC_div4_2,
    const idx_t IC_div4_3, const idx_t IC_div4_4, const idx_t IC_div4_5,
    const idx_t IC_div4_6, const idx_t IC_div4_7, const idx_t IC_div4_8) {
    const idx_t filt_oc_stride = IC_div4_8 + IC_div4;

    __effcc_ignore_memory_order {
        int32_t *output_base = output;
        vec4xi8_t *filter_ptr_a = (vec4xi8_t *)filter;
        vec4xi8_t *filter_ptr_b = filter_ptr_a + filt_oc_stride;

        for (idx_t oc = 0; oc < OC; oc += OC_unroll_2) {
            int32_t b0 = bias[oc];
            int32_t b1 = bias[oc + 1];

            for (idx_t ic = 0; ic < IC; ic += IC_increment_4) {
                /* 9 filter taps for oc */
                vec4xi8_t a00 = filter_ptr_a[0];
                vec4xi8_t a01 = filter_ptr_a[IC_div4];
                vec4xi8_t a02 = filter_ptr_a[IC_div4_2];
                vec4xi8_t a10 = filter_ptr_a[IC_div4_3];
                vec4xi8_t a11 = filter_ptr_a[IC_div4_4];
                vec4xi8_t a12 = filter_ptr_a[IC_div4_5];
                vec4xi8_t a20 = filter_ptr_a[IC_div4_6];
                vec4xi8_t a21 = filter_ptr_a[IC_div4_7];
                vec4xi8_t a22 = filter_ptr_a[IC_div4_8];

                /* 9 filter taps for oc+1 */
                vec4xi8_t b00 = filter_ptr_b[0];
                vec4xi8_t b01 = filter_ptr_b[IC_div4];
                vec4xi8_t b02 = filter_ptr_b[IC_div4_2];
                vec4xi8_t b10 = filter_ptr_b[IC_div4_3];
                vec4xi8_t b11 = filter_ptr_b[IC_div4_4];
                vec4xi8_t b12 = filter_ptr_b[IC_div4_5];
                vec4xi8_t b20 = filter_ptr_b[IC_div4_6];
                vec4xi8_t b21 = filter_ptr_b[IC_div4_7];
                vec4xi8_t b22 = filter_ptr_b[IC_div4_8];

                filter_ptr_a++;
                filter_ptr_b++;

                vec4xi8_t *input_ptr = (vec4xi8_t *)(input_base + ic);
                int32_t *output_ptr = output_base + oc;
                for (idx_t h = 0; h < OH; ++h) {
                    vec4xi8_t in00 = zp_vec, in10 = zp_vec, in20 = zp_vec;

                    bool not_first_row = (h != 0);
                    bool not_last_row = (h != OH - 1);
                    vec4xi8_t in01 =
                        load_guard(not_first_row, input_ptr[0], zp_vec);
                    vec4xi8_t in11 = input_ptr[IW_IC_div4];
                    vec4xi8_t in21 = load_guard(
                        not_last_row, input_ptr[IW_IC_div4_2], zp_vec);

                    vec4xi8_t *input_ptr_0 = input_ptr + IC_div4;
                    vec4xi8_t *input_ptr_1 = input_ptr_0 + IW_IC_div4;
                    vec4xi8_t *input_ptr_2 = input_ptr_1 + IW_IC_div4;
                    input_ptr += IW_IC_div4;

                    for (idx_t w = 0; w < OW; ++w) {
                        bool not_last_col = (w != OW - 1);
                        vec4xi8_t in02 =
                            load_guard((not_first_row && not_last_col),
                                       *input_ptr_0, zp_vec);
                        vec4xi8_t in12 =
                            load_guard(not_last_col, *input_ptr_1, zp_vec);
                        vec4xi8_t in22 =
                            load_guard((not_last_row && not_last_col),
                                       *input_ptr_2, zp_vec);

                        int32_t s00 = sdot_vec4xi8(in00, a00);
                        int32_t s01 = sdot_vec4xi8(in01, a01);
                        int32_t s02 = sdot_vec4xi8(in02, a02);
                        int32_t s10 = sdot_vec4xi8(in10, a10);
                        int32_t s11 = sdot_vec4xi8(in11, a11);
                        int32_t s12 = sdot_vec4xi8(in12, a12);
                        int32_t s20 = sdot_vec4xi8(in20, a20);
                        int32_t s21 = sdot_vec4xi8(in21, a21);
                        int32_t s22 = sdot_vec4xi8(in22, a22);

                        int32_t t00 = sdot_vec4xi8(in00, b00);
                        int32_t t01 = sdot_vec4xi8(in01, b01);
                        int32_t t02 = sdot_vec4xi8(in02, b02);
                        int32_t t10 = sdot_vec4xi8(in10, b10);
                        int32_t t11 = sdot_vec4xi8(in11, b11);
                        int32_t t12 = sdot_vec4xi8(in12, b12);
                        int32_t t20 = sdot_vec4xi8(in20, b20);
                        int32_t t21 = sdot_vec4xi8(in21, b21);
                        int32_t t22 = sdot_vec4xi8(in22, b22);

                        int32_t sum0 =
                            add_tree_8(s00, s01, s02, s10, s11, s12, s20, s21) +
                            s22;
                        int32_t sum1 =
                            add_tree_8(t00, t01, t02, t10, t11, t12, t20, t21) +
                            t22;

                        int32_t run0 = likely(ic != 0) ? output_ptr[0] : b0;
                        int32_t run1 = likely(ic != 0) ? output_ptr[1] : b1;

                        output_ptr[0] = run0 + sum0;
                        output_ptr[1] = run1 + sum1;

                        in00 = in01;
                        in10 = in11;
                        in20 = in21;

                        in01 = in02;
                        in11 = in12;
                        in21 = in22;

                        input_ptr_0 += IC_div4;
                        input_ptr_1 += IC_div4;
                        input_ptr_2 += IC_div4;

                        output_ptr += OC;
                    }
                }
            }
            filter_ptr_a += IC_div4_8 + filt_oc_stride;
            filter_ptr_b += IC_div4_8 + filt_oc_stride;
        }
    }
}

/*
 * 1D systolic × 4 IC groups × 2 OC = 24 sdots per iteration.
 * Optimized reduction: add_tree_8 + add_tree_4 per OC channel.
 * Requires IC_div4 >= 4 (IC >= 16).
 */
__efficient__ void conv2d_stride_1_pad_1_1d_ic4_oc2(
    const int8_t *input, idx_t IW, const int8_t *filter, idx_t OC, idx_t IC,
    int8_t zp, const int32_t *bias, int32_t *output, idx_t OH, idx_t OW,

    const idx_t IC_div4, const idx_t IW_IC_div4, vec4xi8_t zp_vec,
    const idx_t filt_oc_stride, const idx_t filt_kw_stride) {
    const idx_t IW_IC = IW * IC;
    const idx_t filt_oc_stride_2 = filt_oc_stride + filt_oc_stride;

    __effcc_ignore_memory_order {
        int32_t *output_base = output;
        vec4xi8_t *filt_oc_a = (vec4xi8_t *)filter;
        vec4xi8_t *filt_oc_b = filt_oc_a + filt_oc_stride;

        for (idx_t oc = 0; oc < OC; oc += OC_unroll_2) {
            int32_t b0 = bias[oc];
            int32_t b1 = bias[oc + 1];

            const int8_t *input_row_base = input - IW_IC;

            vec4xi8_t *filt_kh_a = filt_oc_a;
            vec4xi8_t *filt_kh_b = filt_oc_b;

            for (idx_t kh = 0; kh < 3; ++kh) {
                vec4xi8_t *fptr_a = filt_kh_a;
                vec4xi8_t *fptr_b = filt_kh_b;

                for (idx_t ic = 0; ic < IC; ic += IC_increment_16) {
                    /* 12 filter taps for oc: 3 kw × 4 IC groups */
                    vec4xi8_t a0_0 = fptr_a[0];
                    vec4xi8_t a0_1 = fptr_a[1];
                    vec4xi8_t a0_2 = fptr_a[2];
                    vec4xi8_t a0_3 = fptr_a[3];
                    vec4xi8_t a1_0 = fptr_a[IC_div4];
                    vec4xi8_t a1_1 = fptr_a[IC_div4 + 1];
                    vec4xi8_t a1_2 = fptr_a[IC_div4 + 2];
                    vec4xi8_t a1_3 = fptr_a[IC_div4 + 3];
                    vec4xi8_t a2_0 = fptr_a[IC_div4 + IC_div4];
                    vec4xi8_t a2_1 = fptr_a[IC_div4 + IC_div4 + 1];
                    vec4xi8_t a2_2 = fptr_a[IC_div4 + IC_div4 + 2];
                    vec4xi8_t a2_3 = fptr_a[IC_div4 + IC_div4 + 3];

                    /* 12 filter taps for oc+1 */
                    vec4xi8_t b0_0 = fptr_b[0];
                    vec4xi8_t b0_1 = fptr_b[1];
                    vec4xi8_t b0_2 = fptr_b[2];
                    vec4xi8_t b0_3 = fptr_b[3];
                    vec4xi8_t b1_0 = fptr_b[IC_div4];
                    vec4xi8_t b1_1 = fptr_b[IC_div4 + 1];
                    vec4xi8_t b1_2 = fptr_b[IC_div4 + 2];
                    vec4xi8_t b1_3 = fptr_b[IC_div4 + 3];
                    vec4xi8_t b2_0 = fptr_b[IC_div4 + IC_div4];
                    vec4xi8_t b2_1 = fptr_b[IC_div4 + IC_div4 + 1];
                    vec4xi8_t b2_2 = fptr_b[IC_div4 + IC_div4 + 2];
                    vec4xi8_t b2_3 = fptr_b[IC_div4 + IC_div4 + 3];

                    fptr_a += 4;
                    fptr_b += 4;

                    vec4xi8_t *input_ptr = (vec4xi8_t *)(input_row_base + ic);
                    int32_t *output_ptr = output_base + oc;

                    for (idx_t h = 0; h < OH; ++h) {
                        bool row_valid =
                            !((h == 0 && kh == 0) || (kh == 2 && h == OH - 1));

                        /* Init sliding window: left = zp, center loaded */
                        vec4xi8_t in_L0 = zp_vec, in_L1 = zp_vec;
                        vec4xi8_t in_L2 = zp_vec, in_L3 = zp_vec;
                        vec4xi8_t in_C0 =
                            load_guard(row_valid, input_ptr[0], zp_vec);
                        vec4xi8_t in_C1 =
                            load_guard(row_valid, input_ptr[1], zp_vec);
                        vec4xi8_t in_C2 =
                            load_guard(row_valid, input_ptr[2], zp_vec);
                        vec4xi8_t in_C3 =
                            load_guard(row_valid, input_ptr[3], zp_vec);

                        vec4xi8_t *right_ptr = input_ptr + IC_div4;
                        input_ptr += IW_IC_div4;

                        for (idx_t w = 0; w < OW; ++w) {
                            bool valid_right = row_valid && (w != OW - 1);

                            vec4xi8_t in_R0 =
                                load_guard(valid_right, right_ptr[0], zp_vec);
                            vec4xi8_t in_R1 =
                                load_guard(valid_right, right_ptr[1], zp_vec);
                            vec4xi8_t in_R2 =
                                load_guard(valid_right, right_ptr[2], zp_vec);
                            vec4xi8_t in_R3 =
                                load_guard(valid_right, right_ptr[3], zp_vec);

                            /* 24 sdots */
                            int32_t p00 = sdot_vec4xi8(in_L0, a0_0);
                            int32_t p01 = sdot_vec4xi8(in_L1, a0_1);
                            int32_t p02 = sdot_vec4xi8(in_L2, a0_2);
                            int32_t p03 = sdot_vec4xi8(in_L3, a0_3);
                            int32_t p10 = sdot_vec4xi8(in_C0, a1_0);
                            int32_t p11 = sdot_vec4xi8(in_C1, a1_1);
                            int32_t p12 = sdot_vec4xi8(in_C2, a1_2);
                            int32_t p13 = sdot_vec4xi8(in_C3, a1_3);
                            int32_t p20 = sdot_vec4xi8(in_R0, a2_0);
                            int32_t p21 = sdot_vec4xi8(in_R1, a2_1);
                            int32_t p22 = sdot_vec4xi8(in_R2, a2_2);
                            int32_t p23 = sdot_vec4xi8(in_R3, a2_3);

                            int32_t q00 = sdot_vec4xi8(in_L0, b0_0);
                            int32_t q01 = sdot_vec4xi8(in_L1, b0_1);
                            int32_t q02 = sdot_vec4xi8(in_L2, b0_2);
                            int32_t q03 = sdot_vec4xi8(in_L3, b0_3);
                            int32_t q10 = sdot_vec4xi8(in_C0, b1_0);
                            int32_t q11 = sdot_vec4xi8(in_C1, b1_1);
                            int32_t q12 = sdot_vec4xi8(in_C2, b1_2);
                            int32_t q13 = sdot_vec4xi8(in_C3, b1_3);
                            int32_t q20 = sdot_vec4xi8(in_R0, b2_0);
                            int32_t q21 = sdot_vec4xi8(in_R1, b2_1);
                            int32_t q22 = sdot_vec4xi8(in_R2, b2_2);
                            int32_t q23 = sdot_vec4xi8(in_R3, b2_3);

                            /* Reduce: add_tree_8 for first 8, add_tree_4 for
                             * last 4, combine */
                            int32_t sum0 = add_tree_8(p00, p01, p02, p03, p10,
                                                      p11, p12, p13) +
                                           add_tree_4(p20, p21, p22, p23);
                            int32_t sum1 = add_tree_8(q00, q01, q02, q03, q10,
                                                      q11, q12, q13) +
                                           add_tree_4(q20, q21, q22, q23);

                            int32_t run0 =
                                likely(ic != 0 || kh != 0) ? output_ptr[0] : b0;
                            int32_t run1 =
                                likely(ic != 0 || kh != 0) ? output_ptr[1] : b1;

                            output_ptr[0] = run0 + sum0;
                            output_ptr[1] = run1 + sum1;

                            in_L0 = in_C0;
                            in_L1 = in_C1;
                            in_L2 = in_C2;
                            in_L3 = in_C3;
                            in_C0 = in_R0;
                            in_C1 = in_R1;
                            in_C2 = in_R2;
                            in_C3 = in_R3;

                            right_ptr += IC_div4;
                            output_ptr += OC;
                        }
                    }
                }

                input_row_base += IW_IC;
                filt_kh_a += filt_kw_stride;
                filt_kh_b += filt_kw_stride;
            }

            filt_oc_a += filt_oc_stride_2;
            filt_oc_b += filt_oc_stride_2;
        }
    }
}

void conv2d_stride_1_pad_1_impl(
    // ==== Input (i8, NHWC) ====
    const int8_t *input, idx_t IW,

    // ==== Filter (i8, [OC,KH,KW,IC]) ====
    const int8_t *filter, idx_t OC, idx_t IC,

    // ==== Zero-point ====
    int8_t zp,

    // ==== Bias (i32, length = OC) ====
    const int32_t *bias,

    // ==== Output (i8, NHWC) ====
    int32_t *output, idx_t OH, idx_t OW) {
    const idx_t IC_div4 = IC / 4;
    const idx_t IW_IC_div4 = IW * IC_div4;
    const idx_t IW_IC_div4_2 = IW_IC_div4 * 2;
    vec4xi8_t zp_vec = {zp, zp, zp, zp};
    const int8_t *input_base = input - IW * IC;
    const idx_t IC_div4_2 = IC_div4 * 2;
    const idx_t IC_div4_3 = IC_div4 * 3;
    const idx_t IC_div4_4 = IC_div4_2 + IC_div4_2;
    const idx_t IC_div4_5 = IC_div4_2 + IC_div4_3;
    const idx_t IC_div4_6 = IC_div4_3 + IC_div4_3;
    const idx_t IC_div4_7 = IC_div4_3 + IC_div4_4;
    const idx_t IC_div4_8 = IC_div4_4 + IC_div4_4;

    if (IC_div4 >= 4) {
        const idx_t filt_oc_stride = IC_div4_8 + IC_div4;
        conv2d_stride_1_pad_1_1d_ic4_oc2(input, IW, filter, OC, IC, zp, bias,
                                         output, OH, OW, IC_div4, IW_IC_div4,
                                         zp_vec, filt_oc_stride, IC_div4_3);
    } else {
        conv2d_stride_1_pad_1_impl_w_const_oc_unroll2(
            input, IW, filter, OC, IC, zp, bias, output, OH, OW, IC_div4,
            IW_IC_div4, IW_IC_div4_2, zp_vec, input_base, IC_div4_2, IC_div4_3,
            IC_div4_4, IC_div4_5, IC_div4_6, IC_div4_7, IC_div4_8);
    }
}

__efficient__ void conv2d_bias_zp(const int8_t *filter, idx_t filter_stride_oc,
                                  const int32_t *bias, idx_t bias_dim_oc,
                                  int32_t zp, int32_t *bias_w_zp) {
    for (idx_t i = 0; i < bias_dim_oc; ++i) {
        int32_t bias_zp = 0;
        for (idx_t j = 0; j < filter_stride_oc; j++) {
            bias_zp += zp * (int32_t)filter[j];
        }
        filter += filter_stride_oc;
        bias_w_zp[i] = bias[i] + bias_zp;
    }
}

void conv2d_stride_1_pad_1(
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

    // ==== Pad0 (i8, length = 1) ====
    const int8_t *pad0_ptr, idx_t pad0_offset, idx_t pad0_dim_1,
    idx_t pad0_stride_1,

    // ==== Pad1 (i8, length = 1) ====
    const int8_t *pad1_ptr, idx_t pad1_offset, idx_t pad1_dim_1,
    idx_t pad1_stride_1,

    // ==== Output (i8, NHWC) ====
    int32_t *out_base, idx_t out_offset, idx_t out_dim_n, idx_t out_dim_h,
    idx_t out_dim_w, idx_t out_dim_c, idx_t out_stride_n, idx_t out_stride_h,
    idx_t out_stride_w, idx_t out_stride_c) {
    int32_t bias_w_zp[bias_dim_oc];
    conv2d_bias_zp(filter_base + filter_offset, filter_stride_oc,
                   bias_base + bias_offset, bias_dim_oc, -(int32_t)pad0_ptr[0],
                   bias_w_zp);
    conv2d_stride_1_pad_1_impl(in_base + in_offset, in_dim_w,
                               filter_base + filter_offset, filter_dim_oc,
                               filter_dim_ic, pad0_ptr[0], bias_w_zp,
                               out_base + out_offset, out_dim_h, out_dim_w);
}

__efficient__ void conv2d_stride_2_pad_corner_w_const_oc_unroll2(
    const int8_t *input, idx_t IW, const int8_t *filter, idx_t OC, idx_t IC,
    int8_t zp, const int32_t *bias, int32_t *output, idx_t OH, idx_t OW,

    const idx_t IC_div4, const idx_t IW_IC_div4, const idx_t IW_IC_div4_2,
    vec4xi8_t zp_vec, const int8_t *input_base, const idx_t IC_div4_2,
    const idx_t IC_div4_3, const idx_t IC_div4_4, const idx_t IC_div4_5,
    const idx_t IC_div4_6, const idx_t IC_div4_7, const idx_t IC_div4_8,
    const idx_t filt_oc_stride, const idx_t IC_div4_17) {
    __effcc_ignore_memory_order {
        int32_t *output_base = output;
        vec4xi8_t *filter_ptr_a = (vec4xi8_t *)filter;
        vec4xi8_t *filter_ptr_b = filter_ptr_a + filt_oc_stride;

        for (idx_t oc = 0; oc < OC; oc += OC_unroll_2) {
            int32_t b0 = bias[oc];
            int32_t b1 = bias[oc + 1];

            for (idx_t ic = 0; ic < IC; ic += IC_increment_4) {
                /* 9 filter taps for oc */
                vec4xi8_t a00 = filter_ptr_a[0];
                vec4xi8_t a01 = filter_ptr_a[IC_div4];
                vec4xi8_t a02 = filter_ptr_a[IC_div4_2];
                vec4xi8_t a10 = filter_ptr_a[IC_div4_3];
                vec4xi8_t a11 = filter_ptr_a[IC_div4_4];
                vec4xi8_t a12 = filter_ptr_a[IC_div4_5];
                vec4xi8_t a20 = filter_ptr_a[IC_div4_6];
                vec4xi8_t a21 = filter_ptr_a[IC_div4_7];
                vec4xi8_t a22 = filter_ptr_a[IC_div4_8];

                /* 9 filter taps for oc+1 */
                vec4xi8_t b00 = filter_ptr_b[0];
                vec4xi8_t b01 = filter_ptr_b[IC_div4];
                vec4xi8_t b02 = filter_ptr_b[IC_div4_2];
                vec4xi8_t b10 = filter_ptr_b[IC_div4_3];
                vec4xi8_t b11 = filter_ptr_b[IC_div4_4];
                vec4xi8_t b12 = filter_ptr_b[IC_div4_5];
                vec4xi8_t b20 = filter_ptr_b[IC_div4_6];
                vec4xi8_t b21 = filter_ptr_b[IC_div4_7];
                vec4xi8_t b22 = filter_ptr_b[IC_div4_8];

                filter_ptr_a++;
                filter_ptr_b++;

                vec4xi8_t *input_ptr = (vec4xi8_t *)(input_base + ic);
                int32_t *output_ptr = output_base + oc;
                for (idx_t h = 0; h < OH; ++h) {
                    bool not_last_row = (h != OH - 1);

                    vec4xi8_t in00 = input_ptr[0];
                    vec4xi8_t in10 = input_ptr[IW_IC_div4];
                    vec4xi8_t in20 = load_guard(
                        not_last_row, input_ptr[IW_IC_div4_2], zp_vec);

                    vec4xi8_t *input_ptr_0 = input_ptr + IC_div4;
                    vec4xi8_t *input_ptr_1 = input_ptr_0 + IW_IC_div4;
                    vec4xi8_t *input_ptr_2 = input_ptr_1 + IW_IC_div4;
                    input_ptr += IW_IC_div4_2;

                    for (idx_t w = 0; w < OW; ++w) {
                        bool not_last_col = (w != OW - 1);
                        vec4xi8_t in01 = input_ptr_0[0];
                        vec4xi8_t in02 = load_guard(
                            not_last_col, input_ptr_0[IC_div4], zp_vec);
                        vec4xi8_t in11 = input_ptr_1[0];
                        vec4xi8_t in12 = load_guard(
                            not_last_col, input_ptr_1[IC_div4], zp_vec);
                        vec4xi8_t in21 =
                            load_guard(not_last_row, input_ptr_2[0], zp_vec);
                        vec4xi8_t in22 =
                            load_guard(not_last_row && not_last_col,
                                       input_ptr_2[IC_div4], zp_vec);

                        int32_t s00 = sdot_vec4xi8(in00, a00);
                        int32_t s01 = sdot_vec4xi8(in01, a01);
                        int32_t s02 = sdot_vec4xi8(in02, a02);
                        int32_t s10 = sdot_vec4xi8(in10, a10);
                        int32_t s11 = sdot_vec4xi8(in11, a11);
                        int32_t s12 = sdot_vec4xi8(in12, a12);
                        int32_t s20 = sdot_vec4xi8(in20, a20);
                        int32_t s21 = sdot_vec4xi8(in21, a21);
                        int32_t s22 = sdot_vec4xi8(in22, a22);

                        int32_t t00 = sdot_vec4xi8(in00, b00);
                        int32_t t01 = sdot_vec4xi8(in01, b01);
                        int32_t t02 = sdot_vec4xi8(in02, b02);
                        int32_t t10 = sdot_vec4xi8(in10, b10);
                        int32_t t11 = sdot_vec4xi8(in11, b11);
                        int32_t t12 = sdot_vec4xi8(in12, b12);
                        int32_t t20 = sdot_vec4xi8(in20, b20);
                        int32_t t21 = sdot_vec4xi8(in21, b21);
                        int32_t t22 = sdot_vec4xi8(in22, b22);

                        int32_t sum0 =
                            add_tree_8(s00, s01, s02, s10, s11, s12, s20, s21) +
                            s22;
                        int32_t sum1 =
                            add_tree_8(t00, t01, t02, t10, t11, t12, t20, t21) +
                            t22;

                        int32_t run0 = likely(ic != 0) ? output_ptr[0] : b0;
                        int32_t run1 = likely(ic != 0) ? output_ptr[1] : b1;

                        output_ptr[0] = run0 + sum0;
                        output_ptr[1] = run1 + sum1;

                        in00 = in02;
                        in10 = in12;
                        in20 = in22;

                        input_ptr_0 += IC_div4_2;
                        input_ptr_1 += IC_div4_2;
                        input_ptr_2 += IC_div4_2;

                        output_ptr += OC;
                    }
                }
            }
            filter_ptr_a += IC_div4_17;
            filter_ptr_b += IC_div4_17;
        }
    }
}

void conv2d_stride_2_pad_corner_impl(const int8_t *input, idx_t IW,
                                     const int8_t *filter, idx_t OC, idx_t IC,
                                     int8_t zp, const int32_t *bias,
                                     int32_t *output, idx_t OH, idx_t OW) {
    const idx_t IC_div4 = IC / 4;
    const idx_t IW_IC_div4 = IW * IC_div4;
    const idx_t IW_IC_div4_2 = IW_IC_div4 * 2;
    vec4xi8_t zp_vec = {zp, zp, zp, zp};
    const int8_t *input_base = input;
    const idx_t IC_div4_2 = IC_div4 * 2;
    const idx_t IC_div4_3 = IC_div4 * 3;
    const idx_t IC_div4_4 = IC_div4_2 + IC_div4_2;
    const idx_t IC_div4_5 = IC_div4_2 + IC_div4_3;
    const idx_t IC_div4_6 = IC_div4_3 + IC_div4_3;
    const idx_t IC_div4_7 = IC_div4_3 + IC_div4_4;
    const idx_t IC_div4_8 = IC_div4_4 + IC_div4_4;
    const idx_t IC_div4_9 = IC_div4_4 + IC_div4_5;
    const idx_t IC_div4_17 = IC_div4_8 + IC_div4_9;

    conv2d_stride_2_pad_corner_w_const_oc_unroll2(
        input, IW, filter, OC, IC, zp, bias, output, OH, OW, IC_div4,
        IW_IC_div4, IW_IC_div4_2, zp_vec, input_base, IC_div4_2, IC_div4_3,
        IC_div4_4, IC_div4_5, IC_div4_6, IC_div4_7, IC_div4_8, IC_div4_9,
        IC_div4_17);
}

void conv2d_stride_2_pad_1_corner(
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

    // ==== Pad0 (i8, length = 1) ====
    const int8_t *pad0_ptr, idx_t pad0_offset, idx_t pad0_dim_1,
    idx_t pad0_stride_1,

    // ==== Pad1 (i8, length = 1) ====
    const int8_t *pad1_ptr, idx_t pad1_offset, idx_t pad1_dim_1,
    idx_t pad1_stride_1,

    // ==== Output (i8, NHWC) ====
    int32_t *out_base, idx_t out_offset, idx_t out_dim_n, idx_t out_dim_h,
    idx_t out_dim_w, idx_t out_dim_c, idx_t out_stride_n, idx_t out_stride_h,
    idx_t out_stride_w, idx_t out_stride_c) {
    int32_t bias_w_zp[bias_dim_oc];
    conv2d_bias_zp(filter_base + filter_offset, filter_stride_oc,
                   bias_base + bias_offset, bias_dim_oc, -(int32_t)pad0_ptr[0],
                   bias_w_zp);
    conv2d_stride_2_pad_corner_impl(
        in_base + in_offset, in_dim_w, filter_base + filter_offset,
        filter_dim_oc, filter_dim_ic, pad0_ptr[0], bias_w_zp,
        out_base + out_offset, out_dim_h, out_dim_w);
}

__efficient__ void pad_channel_to_div4(const int8_t *input, idx_t N, idx_t H,
                                       idx_t W, idx_t C, vec4xi8_t *output) {
    const int8_t *input_ptr = input;
    int32_t *output_ptr = (int32_t *)output;
    for (idx_t n = 0; n < N; n++) {
        for (idx_t h = 0; h < H; h++) {
            for (idx_t w = 0; w < W; w++) {
                idx_t remainder = C;
                for (idx_t c = 0; (c + CH_increment_4) < C;
                     c += CH_increment_4) {
                    *output_ptr++ = *((int32_t *)input_ptr);
                    input_ptr += CH_increment_4;
                    remainder = C - c;
                }
                int32_t val = 0;
                for (idx_t r = 0; r < remainder; r++) {
                    val |= ((int32_t)*input_ptr & 0xFF) << (r * 8);
                    input_ptr++;
                }
                if (remainder) {
                    *output_ptr++ = val;
                }
            }
        }
    }
}

void conv2d_ch_3_stride_1_pad_1(
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

    // ==== Pad0 (i8, length = 1) ====
    const int8_t *pad0_ptr, idx_t pad0_offset, idx_t pad0_dim_1,
    idx_t pad0_stride_1,

    // ==== Pad1 (i8, length = 1) ====
    const int8_t *pad1_ptr, idx_t pad1_offset, idx_t pad1_dim_1,
    idx_t pad1_stride_1,

    // ==== Output (i8, NHWC) ====
    int32_t *out_base, idx_t out_offset, idx_t out_dim_n, idx_t out_dim_h,
    idx_t out_dim_w, idx_t out_dim_c, idx_t out_stride_n, idx_t out_stride_h,
    idx_t out_stride_w, idx_t out_stride_c) {
    int32_t bias_w_zp[bias_dim_oc];
    conv2d_bias_zp(filter_base + filter_offset, filter_stride_oc,
                   bias_base + bias_offset, bias_dim_oc, -(int32_t)pad0_ptr[0],
                   bias_w_zp);

    vec4xi8_t *padded_input =
        (vec4xi8_t *)malloc(in_dim_n * in_dim_h * in_dim_w * sizeof(vec4xi8_t));
    pad_channel_to_div4(in_base + in_offset, in_dim_n, in_dim_h, in_dim_w,
                        in_dim_c, padded_input);

    vec4xi8_t *padded_filter = (vec4xi8_t *)malloc(
        filter_dim_oc * filter_dim_kh * filter_dim_kw * sizeof(vec4xi8_t));
    pad_channel_to_div4(filter_base + filter_offset, filter_dim_oc,
                        filter_dim_kh, filter_dim_kw, filter_dim_ic,
                        padded_filter);

    conv2d_stride_1_pad_1_impl((int8_t *)padded_input, in_dim_w,
                               (int8_t *)padded_filter, filter_dim_oc,
                               filter_dim_ic + 1, pad0_ptr[0], bias_w_zp,
                               out_base + out_offset, out_dim_h, out_dim_w);

    free(padded_input);
    free(padded_filter);
}

/*
 * 1x1 conv with stride 2, no padding.
 * Simple dot product per output pixel: output[oh,ow,oc] = bias[oc] +
 * sum_ic(in[2*oh,2*ow,ic]*filt[oc,ic]) Uses sdot_vec4xi8 for IC reduction.
 * Processes 2 OC at once.
 */
__efficient__ void conv1x1_stride2_impl(const int8_t *input, idx_t in_stride_h,
                                        idx_t in_stride_w, const int8_t *filter,
                                        idx_t OC, idx_t IC, const int32_t *bias,
                                        int32_t *output, idx_t OH, idx_t OW) {
    const idx_t IC_div4 = IC / 4;
    const idx_t filt_oc_stride = IC_div4;

    __effcc_ignore_memory_order {
        int32_t *output_base = output;
        vec4xi8_t *filt_a = (vec4xi8_t *)filter;
        vec4xi8_t *filt_b = filt_a + filt_oc_stride;

        for (idx_t oc = 0; oc < OC; oc += OC_unroll_2) {
            int32_t b0 = bias[oc];
            int32_t b1 = bias[oc + 1];

            const int8_t *inp_row = input;
            int32_t *out_ptr = output_base + oc;

            for (idx_t h = 0; h < OH; ++h) {
                const int8_t *inp_col = inp_row;

                for (idx_t w = 0; w < OW; ++w) {
                    int32_t acc0 = b0;
                    int32_t acc1 = b1;
                    vec4xi8_t *inp_v = (vec4xi8_t *)inp_col;
                    vec4xi8_t *fa = filt_a;
                    vec4xi8_t *fb = filt_b;

                    for (idx_t ic4 = 0; ic4 < IC_div4; ++ic4) {
                        acc0 += sdot_vec4xi8(*inp_v, *fa);
                        acc1 += sdot_vec4xi8(*inp_v, *fb);
                        inp_v++;
                        fa++;
                        fb++;
                    }

                    out_ptr[0] = acc0;
                    out_ptr[1] = acc1;

                    inp_col += in_stride_w + in_stride_w;
                    out_ptr += OC;
                }

                inp_row += in_stride_h + in_stride_h;
            }

            filt_a += filt_oc_stride + filt_oc_stride;
            filt_b += filt_oc_stride + filt_oc_stride;
        }
    }
}

void conv1x1_stride2(const int8_t *in_base, idx_t in_offset, idx_t in_dim_n,
                     idx_t in_dim_h, idx_t in_dim_w, idx_t in_dim_c,
                     idx_t in_stride_n, idx_t in_stride_h, idx_t in_stride_w,
                     idx_t in_stride_c,

                     const int8_t *filter_base, idx_t filter_offset,
                     idx_t filter_dim_oc, idx_t filter_dim_kh,
                     idx_t filter_dim_kw, idx_t filter_dim_ic,
                     idx_t filter_stride_oc, idx_t filter_stride_h,
                     idx_t filter_stride_w, idx_t filter_stride_ic,

                     const int32_t *bias_base, idx_t bias_offset,
                     idx_t bias_dim_oc, idx_t bias_stride_1,

                     const int8_t *pad0_ptr, idx_t pad0_offset,
                     idx_t pad0_dim_1, idx_t pad0_stride_1,

                     const int8_t *pad1_ptr, idx_t pad1_offset,
                     idx_t pad1_dim_1, idx_t pad1_stride_1,

                     int32_t *out_base, idx_t out_offset, idx_t out_dim_n,
                     idx_t out_dim_h, idx_t out_dim_w, idx_t out_dim_c,
                     idx_t out_stride_n, idx_t out_stride_h, idx_t out_stride_w,
                     idx_t out_stride_c) {
    /* Fold zero-point into bias: bias_adj = bias - zp * sum(filter_row) */
    int32_t bias_w_zp[filter_dim_oc];
    conv2d_bias_zp(filter_base + filter_offset, filter_stride_oc,
                   bias_base + bias_offset, filter_dim_oc,
                   -(int32_t)pad0_ptr[0], bias_w_zp);

    /* Use in_stride_h/w (not in_dim_w*IC) because input may be a slice */
    conv1x1_stride2_impl(in_base + in_offset, in_stride_h, in_stride_w,
                         filter_base + filter_offset, filter_dim_oc,
                         filter_dim_ic, bias_w_zp, out_base + out_offset,
                         out_dim_h, out_dim_w);
}
