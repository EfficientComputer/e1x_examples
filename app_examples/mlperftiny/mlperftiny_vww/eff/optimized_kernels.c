// <<AUTOBENCH>> efficient_e1x
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include "effintrinsic.h"

#define likely(x) __builtin_expect(!!(x), 1)

#define LOAD_FILTER_3x3(TYPE, filter_ptr, STRIDE) \
    TYPE f00 = (filter_ptr)[0 * (STRIDE)];        \
    TYPE f01 = (filter_ptr)[1 * (STRIDE)];        \
    TYPE f02 = (filter_ptr)[2 * (STRIDE)];        \
    TYPE f10 = (filter_ptr)[3 * (STRIDE)];        \
    TYPE f11 = (filter_ptr)[4 * (STRIDE)];        \
    TYPE f12 = (filter_ptr)[5 * (STRIDE)];        \
    TYPE f20 = (filter_ptr)[6 * (STRIDE)];        \
    TYPE f21 = (filter_ptr)[7 * (STRIDE)];        \
    TYPE f22 = (filter_ptr)[8 * (STRIDE)];

// element-wise product over 3x3 hw-spatial grid
#define MUL_3x3()            \
    int32_t s0 = f00 * in00; \
    int32_t s1 = f01 * in01; \
    int32_t s2 = f02 * in02; \
    int32_t s3 = f10 * in10; \
    int32_t s4 = f11 * in11; \
    int32_t s5 = f12 * in12; \
    int32_t s6 = f20 * in20; \
    int32_t s7 = f21 * in21; \
    int32_t s8 = f22 * in22;

#define SDOT_3x3()                         \
    int32_t i00 = sdot_vec4xi8(in00, f00); \
    int32_t i01 = sdot_vec4xi8(in01, f01); \
    int32_t i02 = sdot_vec4xi8(in02, f02); \
    int32_t i10 = sdot_vec4xi8(in10, f10); \
    int32_t i11 = sdot_vec4xi8(in11, f11); \
    int32_t i12 = sdot_vec4xi8(in12, f12); \
    int32_t i20 = sdot_vec4xi8(in20, f20); \
    int32_t i21 = sdot_vec4xi8(in21, f21); \
    int32_t i22 = sdot_vec4xi8(in22, f22);

typedef uint32_t idx_t;
typedef int8_t vec4xi8_t __attribute__((ext_vector_type(4)));
extern int32_t sdot_vec4xi8(vec4xi8_t a, vec4xi8_t b);
extern int32_t add_tree_8(int32_t a, int32_t b, int32_t c, int32_t d, int32_t e,
                          int32_t f, int32_t g, int32_t h);
extern int32_t add_tree_10(int32_t a, int32_t b, int32_t c, int32_t d,
                           int32_t e, int32_t f, int32_t g, int32_t h,
                           int32_t i, int32_t j);
extern vec4xi8_t load_guard(bool cond, vec4xi8_t ld, vec4xi8_t alt);
extern int8_t load_guard_i8(bool cond, int8_t ld, int8_t alt);
vec4xi8_t vec_zero = {0, 0, 0, 0};

__efficient__ void conv2d_stride_2_impl(
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
    __effcc_ignore_memory_order {
        const idx_t IC_div4 = IC >> 2;
        const idx_t IW_IC_div4 = IW * IC_div4;
        vec4xi8_t zp_vec = {zp, zp, zp, zp};
        // const int8_t *input_base = input - IW * IC;
        const int8_t *input_base = input;
        const idx_t IC_div4_2 = IC_div4 << 1;
        const idx_t IW_IC_div4_2 = IW_IC_div4 << 1;

        int32_t *output_base = output;
        vec4xi8_t *filter_ptr = (vec4xi8_t *)filter;
        for (idx_t oc = 0; oc < OC; ++oc) {
            int32_t b = bias[oc];

            for (idx_t ic = 0; ic < IC; ic += 4) {
                LOAD_FILTER_3x3(vec4xi8_t, filter_ptr, IC_div4);
                filter_ptr++;

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
                        vec4xi8_t in11 = input_ptr_1[0];
                        vec4xi8_t in21 =
                            load_guard(not_last_row, input_ptr_2[0], zp_vec);

                        vec4xi8_t in02 = load_guard(
                            not_last_col, input_ptr_0[IC_div4], zp_vec);
                        vec4xi8_t in12 = load_guard(
                            not_last_col, input_ptr_1[IC_div4], zp_vec);
                        vec4xi8_t in22 =
                            load_guard(not_last_row && not_last_col,
                                       input_ptr_2[IC_div4], zp_vec);

                        SDOT_3x3();

                        int32_t sum =
                            add_tree_8(i00, i01, i02, i10, i11, i12, i20, i21) +
                            i22;
                        int32_t running = likely(ic != 0) ? *output_ptr : b;

                        *output_ptr = running + sum;

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
            filter_ptr += 8 * IC_div4;
        }
    }
}

__efficient__ void pad_channel_to_div4(const int8_t *input, idx_t N, idx_t H,
                                       idx_t W, idx_t C, vec4xi8_t *output) {
    __effcc_ignore_memory_order {
        const int8_t *input_ptr = input;
        int32_t *output_ptr = (int32_t *)output;
        for (idx_t n = 0; n < N; n++) {
            for (idx_t h = 0; h < H; h++) {
                for (idx_t w = 0; w < W; w++) {
                    idx_t remainder = C;
                    for (idx_t c = 0; (c + 4) < C; c += 4) {
                        *output_ptr++ = *((int32_t *)input_ptr);
                        input_ptr += 4;
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
}

__effcc_rip_exact void conv2d_dw_stride_1_core(
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
    __effcc_ignore_memory_order {
        const idx_t IC_div4 = IC;
        const idx_t IW_IC = IW * IC;
        const idx_t IW_IC_2 = (IW_IC << 1);
        const int8_t *input_base = input - IW * IC;

        int32_t *output_base = output;
        const int8_t *filter_ptr = filter;

        for (idx_t ic = 0; ic < IC; ic++) {
            LOAD_FILTER_3x3(int8_t, filter_ptr, IC);
            filter_ptr++;

            int8_t *input_ptr = (int8_t *)(input_base + ic);
            int32_t *output_ptr = output_base + ic;
            int32_t bias_oc = bias[ic];

            for (idx_t h = 0; h < OH; ++h) {
                int8_t in00 = zp, in10 = zp, in20 = zp;

                bool not_first_row = (h != 0);
                bool not_last_row = (h != OH - 1);

                int8_t in01 = load_guard_i8((not_first_row), input_ptr[0], zp);
                int8_t in11 = input_ptr[IW_IC];
                int8_t in21 =
                    load_guard_i8((not_last_row), input_ptr[IW_IC_2], zp);

                int8_t *input_ptr_0 = input_ptr + IC;
                int8_t *input_ptr_1 = input_ptr_0 + IW_IC;
                int8_t *input_ptr_2 = input_ptr_1 + IW_IC;
                input_ptr += IW_IC;

                for (idx_t w = 0; w < OW; w++) {
                    bool not_last_col = (w != OW - 1);

                    int8_t in02 = load_guard_i8((not_first_row && not_last_col),
                                                *input_ptr_0, zp);
                    int8_t in12 =
                        load_guard_i8((not_last_col), *input_ptr_1, zp);
                    int8_t in22 = load_guard_i8((not_last_row && not_last_col),
                                                *input_ptr_2, zp);

                    // element-wise product  s = in * f
                    MUL_3x3();

                    *output_ptr = add_tree_10(s0, s1, s2, s3, s4, s5, s6, s7,
                                              s8, bias_oc);

                    in00 = in01;
                    in10 = in11;
                    in20 = in21;

                    in01 = in02;
                    in11 = in12;
                    in21 = in22;

                    input_ptr_0 += IC;
                    input_ptr_1 += IC;
                    input_ptr_2 += IC;

                    output_ptr += OC;
                }
            }
        }
    }
}

__effcc_rip_exact void conv2d_dw_stride_2_core(
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
    __effcc_ignore_memory_order {
        const idx_t IW_IC = IW * IC;
        const idx_t IW_IC_2 = (IW_IC << 1);
        const int8_t *input_base = input;  // No top-padding for stride-2.
                                           // So do not adjust the input base
        const idx_t IC_2 = IC << 1;

        int32_t *output_base = output;
        const int8_t *filter_ptr = filter;

        for (idx_t ic = 0; ic < IC; ic++) {
            LOAD_FILTER_3x3(int8_t, filter_ptr, IC);
            filter_ptr++;

            int8_t *input_ptr = (int8_t *)(input_base + ic);
            int32_t *output_ptr = output_base + ic;
            int32_t bias_oc = bias[ic];

            for (idx_t h = 0; h < OH; ++h) {
                bool not_last_row = (h != OH - 1);

                int8_t in00 = input_ptr[0];
                int8_t in10 = input_ptr[IW_IC];
                int8_t in20 =
                    load_guard_i8((not_last_row), input_ptr[IW_IC_2], zp);

                int8_t *input_ptr_0 = input_ptr + IC;
                int8_t *input_ptr_1 = input_ptr_0 + IW_IC;
                int8_t *input_ptr_2 = input_ptr_1 + IW_IC;
                input_ptr += IW_IC_2;

                for (idx_t w = 0; w < OW; w++) {
                    bool not_last_col = (w != OW - 1);

                    int8_t in01 = input_ptr_0[0];
                    int8_t in02 =
                        load_guard_i8(not_last_col, input_ptr_0[IC], zp);
                    int8_t in11 = input_ptr_1[0];
                    int8_t in12 =
                        load_guard_i8(not_last_col, input_ptr_1[IC], zp);
                    int8_t in21 =
                        load_guard_i8(not_last_row, input_ptr_2[0], zp);
                    int8_t in22 = load_guard_i8(not_last_row && not_last_col,
                                                input_ptr_2[IC], zp);

                    // element-wise product s = in * f
                    MUL_3x3();

                    *output_ptr = add_tree_10(s0, s1, s2, s3, s4, s5, s6, s7,
                                              s8, bias_oc);

                    in00 = in02;
                    in10 = in12;
                    in20 = in22;

                    input_ptr_0 += IC_2;
                    input_ptr_1 += IC_2;
                    input_ptr_2 += IC_2;

                    output_ptr += OC;
                }
            }
        }
    }
}

__efficient__ void conv2d_bias_zp_1(const int8_t *filter,
                                    idx_t filter_stride_oc, const int32_t *bias,
                                    idx_t bias_dim_oc, int32_t zp,
                                    int32_t *bias_w_zp) {
    __effcc_ignore_memory_order {
        for (idx_t i = 0; i < bias_dim_oc; ++i) {
            int32_t bias_zp = 0;
            __effcc_parallel(16) for (idx_t j = 0; j < filter_stride_oc; j++) {
                bias_zp += zp * (int32_t)filter[j];
            }
            filter += filter_stride_oc;
            bias_w_zp[i] = bias[i] + bias_zp;
        }
    }
}

__efficient__ void conv2d_bias_zp_2(const int8_t *filter,
                                    idx_t filter_stride_oc, const int32_t *bias,
                                    idx_t bias_dim_oc, int32_t zp,
                                    int32_t *bias_w_zp) {
    // Assumption for depthwise TOSA-style filter layout:
    // filter is [KH, KW, C, 1] (or [3,3,C,1]) and bias_dim_oc == C.
    // Then the 9 taps for channel i are at:
    //   filter[(kh*KW + kw)*C + i]
    //
    // For 3x3: filter_stride_oc should be 9 * bias_dim_oc (i.e., 9*C).
    __effcc_ignore_memory_order {
        const idx_t KH = 3;
        const idx_t KW = 3;
        const idx_t C = bias_dim_oc;

        for (idx_t i = 0; i < bias_dim_oc; ++i) {
            int32_t sumW = 0;

            // Sum the 3x3 weights for channel i (not contiguous in memory!)
            for (idx_t kh = 0; kh < KH; ++kh) {
                for (idx_t kw = 0; kw < KW; ++kw) {
                    const idx_t idx = (kh * KW + kw) * C + i;  // depthwise M=1
                    sumW += (int32_t)filter[idx];
                }
            }

            // Correct sign: bias' = bias - zp * sumW
            bias_w_zp[i] = bias[i] + zp * sumW;
        }
    }
}

void conv2d_dw_stride_1(
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
    conv2d_bias_zp_2(filter_base + filter_offset, filter_stride_oc,
                     bias_base + bias_offset, bias_dim_oc,
                     -(int32_t)pad0_ptr[0], bias_w_zp);

    conv2d_dw_stride_1_core(in_base + in_offset, in_dim_w,
                            filter_base + filter_offset, filter_dim_kw,
                            filter_dim_kw, pad0_ptr[0], bias_w_zp,
                            out_base + out_offset, out_dim_h, out_dim_w);
}

void conv2d_dw_stride_2(
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
    conv2d_bias_zp_2(filter_base + filter_offset, filter_stride_oc,
                     bias_base + bias_offset, bias_dim_oc,
                     -(int32_t)pad0_ptr[0], bias_w_zp);

    conv2d_dw_stride_2_core(in_base + in_offset, in_dim_w,
                            filter_base + filter_offset, filter_dim_kw,
                            filter_dim_kw, pad0_ptr[0], bias_w_zp,
                            out_base + out_offset, out_dim_h, out_dim_w);
}

void conv2d_ch_3_stride_2(
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
    conv2d_bias_zp_1(filter_base + filter_offset, filter_stride_oc,
                     bias_base + bias_offset, bias_dim_oc,
                     -(int32_t)pad0_ptr[0], bias_w_zp);

    vec4xi8_t *padded_input =
        (vec4xi8_t *)malloc(in_dim_n * in_dim_h * in_dim_w * sizeof(vec4xi8_t));
    pad_channel_to_div4(in_base + in_offset, in_dim_n, in_dim_h, in_dim_w,
                        in_dim_c, padded_input);

    vec4xi8_t *padded_filter = (vec4xi8_t *)malloc(
        filter_dim_oc * filter_dim_kh * filter_dim_kw * sizeof(vec4xi8_t));
    pad_channel_to_div4(filter_base + filter_offset, filter_dim_oc,
                        filter_dim_kh, filter_dim_kw, filter_dim_ic,
                        padded_filter);

    conv2d_stride_2_impl((int8_t *)padded_input, in_dim_w,
                         (int8_t *)padded_filter, filter_dim_oc,
                         filter_dim_ic + 1, pad0_ptr[0], bias_w_zp,
                         out_base + out_offset, out_dim_h, out_dim_w);

    free(padded_input);
    free(padded_filter);
}
