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

#define LOAD_VEC8(PREFIX, PTR)            \
    const vec4xi8_t PREFIX##0 = (PTR)[0]; \
    const vec4xi8_t PREFIX##1 = (PTR)[1]; \
    const vec4xi8_t PREFIX##2 = (PTR)[2]; \
    const vec4xi8_t PREFIX##3 = (PTR)[3]; \
    const vec4xi8_t PREFIX##4 = (PTR)[4]; \
    const vec4xi8_t PREFIX##5 = (PTR)[5]; \
    const vec4xi8_t PREFIX##6 = (PTR)[6]; \
    const vec4xi8_t PREFIX##7 = (PTR)[7];

#define SDOT_VEC8(PREFIX, X, K)                   \
    int32_t PREFIX##0 = sdot_vec4xi8(X##0, K##0); \
    int32_t PREFIX##1 = sdot_vec4xi8(X##1, K##1); \
    int32_t PREFIX##2 = sdot_vec4xi8(X##2, K##2); \
    int32_t PREFIX##3 = sdot_vec4xi8(X##3, K##3); \
    int32_t PREFIX##4 = sdot_vec4xi8(X##4, K##4); \
    int32_t PREFIX##5 = sdot_vec4xi8(X##5, K##5); \
    int32_t PREFIX##6 = sdot_vec4xi8(X##6, K##6); \
    int32_t PREFIX##7 = sdot_vec4xi8(X##7, K##7);

typedef uint32_t idx_t;
typedef int8_t vec4xi8_t __attribute__((ext_vector_type(4)));
extern int32_t sdot_vec4xi8(vec4xi8_t a, vec4xi8_t b);
extern int32_t add_tree_8(int32_t a, int32_t b, int32_t c, int32_t d, int32_t e,
                          int32_t f, int32_t g, int32_t h);
extern vec4xi8_t load_guard(bool cond, vec4xi8_t ld, vec4xi8_t alt);

__efficient__ void conv2d_stride_1_pad_1_impl(
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

    __effcc_ignore_memory_order {
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

                        SDOT_3x3();

                        int32_t sum =
                            add_tree_8(i00, i01, i02, i10, i11, i12, i20, i21) +
                            i22;
                        int32_t running = likely(ic != 0) ? *output_ptr : b;

                        *output_ptr = running + sum;

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
            filter_ptr += 8 * IC_div4;
        }
    }
}

__effcc_rip_exact void conv2d_stride_1_pad_1_impl_32_core(
    const int8_t *inBase, idx_t W, idx_t H, const int8_t *filterBase, idx_t OC,
    idx_t IC,  // IC ignored (assumed 32)
    int8_t zp, const int32_t *biasBase, int32_t *outBase, idx_t OH, idx_t OW) {
    __effcc_ignore_memory_order {
        const idx_t OH_core = OH - 2;
        const idx_t OW_core = OW - 2;

        const idx_t filtAdvance = 9 * IC;
        const idx_t outRowStride = OW * OC;

        for (idx_t oc = 0; oc < OC; ++oc) {
            const int32_t bias_val = biasBase[oc];

            // Filter base for this oc
            const int8_t *filtOcBase = filterBase + oc * filtAdvance;

            // out pointer for h=0 (i.e., output row 1), output col 1, channel
            // oc
            int32_t *outRowPtr = outBase + outRowStride + OC + oc;

            for (idx_t h = 0; h < OH_core; ++h) {
                idx_t kh = 0, kw = 0;
                idx_t filterHWAdvance = 0;

                // ---- accumulate 9 taps directly into output ----
                for (idx_t f = 0; f < 9; ++f) {
                    const int8_t *kTapPtr = filtOcBase + filterHWAdvance;

                    vec4xi8_t *kPos = (vec4xi8_t *)kTapPtr;
                    LOAD_VEC8(k, kPos);  // k0..k7

                    const idx_t in_h = h + kh;
                    const idx_t inHWOffset = (in_h * W + kw) << 5;  // *32
                    const int8_t *inPtr = inBase + inHWOffset;

                    // output accumulator pointer for w sweep
                    int32_t *o = outRowPtr;

                    for (idx_t w = 0; w < OW_core; ++w) {
                        vec4xi8_t *inPos = (vec4xi8_t *)inPtr;
                        LOAD_VEC8(x, inPos);  // x0..x7

                        SDOT_VEC8(s, x, k);  // s0..s7

                        int32_t acc = likely(f != 0) ? *o : bias_val;
                        int32_t sum =
                            add_tree_8(s0, s1, s2, s3, s4, s5, s6, s7);

                        *o = sum + acc;

                        inPtr += 32;  // next pixel (IC=32)
                        o += OC;      // next output pixel (same oc)
                    }

                    // advance (kh,kw) and filter tap
                    kw = (kw == 2) ? 0 : (kw + 1);
                    kh = (kw == 0) ? (kh + 1) : kh;
                    filterHWAdvance += 32;
                }

                // bump to next output row (still shifted by +1 col, same oc)
                outRowPtr += outRowStride;
            }
        }
    }
}

__effcc_rip_exact void conv2d_stride_2_pad_1_impl_32_core(
    // ==== Input (i8, NHWC) ====
    const int8_t *inBase, idx_t W, idx_t H,

    // ==== Filter (i8, [OC,KH,KW,IC]) ====
    const int8_t *filterBase, idx_t OC, idx_t IC,

    // ==== Zero-point ====
    int8_t zp,

    // ==== Bias (i32, length = OC) ====
    const int32_t *biasBase,

    // ==== Output (i8, NHWC) ====
    int32_t *outBase, idx_t OH, idx_t OW) {
    __effcc_ignore_memory_order {
        idx_t filtOcOffset = 0;
        const idx_t filtAdvance = 9 * IC;
        const idx_t outRowStride = OW * OC;

        for (idx_t oc = 0; oc < OC; ++oc) {
            const int32_t bias_val = biasBase[oc];

            // output pointer for start of row h=0 at this oc
            int32_t *outPtrRow = outBase + oc;

            for (idx_t h = 0; h < OH - 1; ++h) {
                int32_t *outPtr = outPtrRow;  // w=0 for this row

                for (idx_t w = 0; w < OW - 1; ++w) {
                    int32_t partial_sum = bias_val;

                    // 3×3 spatial window, flattened
                    idx_t kh = 0, kw = 0;
                    idx_t filterHWAdvance = 0;

                    for (idx_t f = 0; f < 9; ++f) {
                        const idx_t in_h = (h << 1) + kh;
                        const idx_t in_w = (w << 1) + kw;

                        // NOTE: still uses mults for input addressing; only
                        // output indexing removed here
                        const idx_t inHWOffset = (in_h * W + in_w) * IC;
                        const idx_t filtHWOffset =
                            filtOcOffset + filterHWAdvance;

                        vec4xi8_t *inPos = (vec4xi8_t *)(inBase + inHWOffset);
                        vec4xi8_t *kPos =
                            (vec4xi8_t *)(filterBase + filtHWOffset);

                        LOAD_VEC8(x, inPos);  // x0..x7
                        LOAD_VEC8(k, kPos);   // k0..k7
                        SDOT_VEC8(s, x, k);   // s0..s7

                        partial_sum +=
                            add_tree_8(s0, s1, s2, s3, s4, s5, s6, s7);

                        kw = (kw == 2) ? 0 : (kw + 1);
                        kh = (kw == 0) ? (kh + 1) : kh;

                        filterHWAdvance += IC;
                    }

                    *outPtr = partial_sum;  // store (h,w,oc)
                    outPtr += OC;           // bump to next w
                }

                outPtrRow += outRowStride;  // bump to next h
            }

            filtOcOffset += filtAdvance;  // next OC filter set
        }
    }
}

__effcc_rip_exact void conv2d_stride_1_pad_1_impl_32_top(
    // ==== Input (i8, NHWC) ====
    const int8_t *inBase, idx_t W, idx_t H,

    // ==== Filter (i8, [OC,KH,KW,IC]) ====
    const int8_t *filterBase, idx_t OC, idx_t IC,

    // ==== Zero-point ====
    int8_t zp,

    // ==== Bias (i32, length = OC) ====
    const int32_t *biasBase,

    // ==== Output (i8, NHWC) ====
    int32_t *outBase, idx_t OH, idx_t OW) {
    __effcc_ignore_memory_order {
        idx_t filtOcOffset = 0;
        const idx_t filtStartOffset =
            3 * IC;  // skip first row of filter taps (kh=0)
        const idx_t filtAdvance = 9 * IC;

        for (idx_t oc = 0; oc < OC; ++oc) {
            const int32_t bias_val = biasBase[oc];

            // h == 0 always in this kernel
            // first output is at (h=0, w=1, oc)
            int32_t *outPtr = outBase + OC + oc;

            for (idx_t w = 0; w < OW - 2; ++w) {
                int32_t partial_sum = bias_val;

                idx_t kh = 0, kw = 0;
                idx_t filterHWAdvance = 0;

                for (idx_t f = 0; f < 6; ++f) {
                    const idx_t in_h = /*h*/ 0 + kh;
                    const idx_t in_w = w + kw;

                    const idx_t inHWOffset = (in_h * W + in_w) * IC;
                    const idx_t filtHWOffset =
                        filtOcOffset + filtStartOffset + filterHWAdvance;

                    vec4xi8_t *inPos = (vec4xi8_t *)(inBase + inHWOffset);
                    vec4xi8_t *kPos = (vec4xi8_t *)(filterBase + filtHWOffset);

                    LOAD_VEC8(x, inPos);  // x0..x7
                    LOAD_VEC8(k, kPos);   // k0..k7
                    SDOT_VEC8(s, x, k);   // s0..s7

                    partial_sum += add_tree_8(s0, s1, s2, s3, s4, s5, s6, s7);

                    kw = (kw == 2) ? 0 : (kw + 1);
                    kh = (kw == 0) ? (kh + 1) : kh;

                    filterHWAdvance += IC;
                }

                *outPtr = partial_sum;  // store at (0, w+1, oc)
                outPtr += OC;           // bump to next w (still +1 col shift)
            }

            filtOcOffset += filtAdvance;
        }
    }
}

__effcc_rip_exact void conv2d_stride_1_pad_1_impl_32_bottom(
    // ==== Input (i8, NHWC) ====
    const int8_t *inBase, idx_t W, idx_t H,

    // ==== Filter (i8, [OC,KH,KW,IC]) ====
    const int8_t *filterBase, idx_t OC, idx_t IC,

    // ==== Zero-point ====
    int8_t zp,

    // ==== Bias (i32, length = OC) ====
    const int32_t *biasBase,

    // ==== Output (i8, NHWC) ====
    int32_t *outBase, idx_t OH, idx_t OW) {
    __effcc_ignore_memory_order {
        idx_t filtOcOffset = 0;
        const idx_t filtStartOffset = 0;
        const idx_t filtAdvance = 9 * IC;

        // fixed input row for this bottom kernel
        const idx_t h = H - 2;

        // output row we write to (shifted down by 1)
        const idx_t out_h = h + 1;

        // precompute output row stride (elements)
        const idx_t outRowStride = OW * OC;

        for (idx_t oc = 0; oc < OC; ++oc) {
            const int32_t bias_val = biasBase[oc];

            // first output is at (out_h, w=1, oc)
            int32_t *outPtr = outBase + (out_h * outRowStride) + OC + oc;

            for (idx_t w = 0; w < OW - 2; ++w) {
                int32_t partial_sum = bias_val;

                idx_t kh = 0, kw = 0;
                idx_t filterHWAdvance = 0;

                for (idx_t f = 0; f < 6; ++f) {
                    const idx_t in_h = h + kh;
                    const idx_t in_w = w + kw;

                    const idx_t inHWOffset = (in_h * W + in_w) * IC;
                    const idx_t filtHWOffset =
                        filtOcOffset + filtStartOffset + filterHWAdvance;

                    vec4xi8_t *inPos = (vec4xi8_t *)(inBase + inHWOffset);
                    vec4xi8_t *kPos = (vec4xi8_t *)(filterBase + filtHWOffset);

                    LOAD_VEC8(x, inPos);  // x0..x7
                    LOAD_VEC8(k, kPos);   // k0..k7
                    SDOT_VEC8(s, x, k);   // s0..s7

                    partial_sum += add_tree_8(s0, s1, s2, s3, s4, s5, s6, s7);

                    kw = (kw == 2) ? 0 : (kw + 1);
                    kh = (kw == 0) ? (kh + 1) : kh;

                    filterHWAdvance += IC;
                }

                *outPtr = partial_sum;  // store at (out_h, w+1, oc)
                outPtr += OC;           // bump to next w
            }

            filtOcOffset += filtAdvance;
        }
    }
}

__effcc_rip_exact void conv2d_stride_2_pad_1_impl_32_bottom(
    // ==== Input (i8, NHWC) ====
    const int8_t *inBase, idx_t W, idx_t H,

    // ==== Filter (i8, [OC,KH,KW,IC]) ====
    const int8_t *filterBase, idx_t OC, idx_t IC,

    // ==== Zero-point ====
    int8_t zp,

    // ==== Bias (i32, length = OC) ====
    const int32_t *biasBase,

    // ==== Output (i8, NHWC) ====
    int32_t *outBase, idx_t OH, idx_t OW) {
    __effcc_ignore_memory_order {
        idx_t filtOcOffset = 0;
        const idx_t filtStartOffset = 0;
        const idx_t filtAdvance = 9 * IC;

        const idx_t h = H - 2;    // fixed input base row for this kernel
        const idx_t oh = OH - 1;  // fixed output row for this kernel

        const idx_t outRowStride = OW * OC;

        for (idx_t oc = 0; oc < OC; ++oc) {
            const int32_t bias_val = biasBase[oc];

            // first output is at (oh, ow=0, oc)
            int32_t *outPtr = outBase + (oh * outRowStride) + oc;

            for (idx_t ow = 0; ow < OW - 1; ++ow) {
                int32_t partial_sum = bias_val;

                const idx_t w = ow << 1;

                idx_t kh = 0, kw = 0;
                idx_t filterHWAdvance = 0;

                for (idx_t f = 0; f < 6; ++f) {
                    const idx_t in_h = h + kh;
                    const idx_t in_w = w + kw;

                    const idx_t inHWOffset = (in_h * W + in_w) * IC;
                    const idx_t filtHWOffset =
                        filtOcOffset + filtStartOffset + filterHWAdvance;

                    vec4xi8_t *inPos = (vec4xi8_t *)(inBase + inHWOffset);
                    vec4xi8_t *kPos = (vec4xi8_t *)(filterBase + filtHWOffset);

                    LOAD_VEC8(x, inPos);  // x0..x7
                    LOAD_VEC8(k, kPos);   // k0..k7
                    SDOT_VEC8(s, x, k);   // s0..s7

                    partial_sum += add_tree_8(s0, s1, s2, s3, s4, s5, s6, s7);

                    kw = (kw == 2) ? 0 : (kw + 1);
                    kh = (kw == 0) ? (kh + 1) : kh;

                    filterHWAdvance += IC;
                }

                *outPtr = partial_sum;  // store at (oh, ow, oc)
                outPtr += OC;           // bump to next ow
            }

            filtOcOffset += filtAdvance;
        }
    }
}

__effcc_rip_exact void conv2d_stride_1_pad_1_impl_32_left(
    // ==== Input (i8, NHWC) ====
    const int8_t *inBase, idx_t W, idx_t H,

    // ==== Filter (i8, [OC,KH,KW,IC]) ====
    const int8_t *filterBase, idx_t OC, idx_t IC,

    // ==== Zero-point ====
    int8_t zp,

    // ==== Bias (i32, length = OC) ====
    const int32_t *biasBase,

    // ==== Output (i8, NHWC) ====
    int32_t *outBase, idx_t OH, idx_t OW) {
    __effcc_ignore_memory_order {
        idx_t filtOcOffset = 0;
        idx_t filtAdvance = 9 * IC;

        for (idx_t oc = 0; oc < OC; ++oc) {
            int32_t bias_val = biasBase[oc];

            idx_t w = 0;
            for (idx_t h = 0; h < OH - 2; ++h) {
                int32_t partial_sum = bias_val;

                // 3×3 spatial window, flattened
                idx_t kh = 0, kw = 0;
                idx_t filterHWAdvance = 0;

                for (idx_t f = 0; f < 6; f++) {
                    idx_t in_h = h + kh;
                    idx_t in_w = w + kw;

                    // Base input index for (n=0, in_h, in_w, 0)
                    idx_t inHWOffset = (in_h * W + in_w) * IC;
                    // Base filter index for (oc, kh, kw, 0)
                    idx_t filtHWOffset =
                        filtOcOffset + ((kh * 3) + kw + 1) * IC;

                    vec4xi8_t *inPos = (vec4xi8_t *)(inBase + inHWOffset);
                    vec4xi8_t *kPos = (vec4xi8_t *)(filterBase + filtHWOffset);

                    LOAD_VEC8(x, inPos);  // x0..x7

                    LOAD_VEC8(k, kPos);  // k0..k7

                    SDOT_VEC8(s, x, k);  // s0..s7

                    partial_sum += add_tree_8(s0, s1, s2, s3, s4, s5, s6, s7);

                    kw = (kw == 1) ? 0 : kw + 1;
                    kh = (kw == 0) ? kh + 1 : kh;

                    // filterHWAdvance += IC;
                }

                // Store output at (n=0, h, w, oc)
                idx_t outIndex = ((h + 1) * OW) * OC + oc;
                outBase[outIndex] = partial_sum;
            }

            filtOcOffset += filtAdvance;  // next OC filter set
        }
    }
}

__effcc_rip_exact void conv2d_stride_1_pad_1_impl_32_right(
    // ==== Input (i8, NHWC) ====
    const int8_t *inBase, idx_t W, idx_t H,

    // ==== Filter (i8, [OC,KH,KW,IC]) ====
    const int8_t *filterBase, idx_t OC, idx_t IC,

    // ==== Zero-point ====
    int8_t zp,

    // ==== Bias (i32, length = OC) ====
    const int32_t *biasBase,

    // ==== Output (i8, NHWC) ====
    int32_t *outBase, idx_t OH, idx_t OW) {
    __effcc_ignore_memory_order {
        idx_t filtOcOffset = 0;
        const idx_t filtAdvance = 9 * IC;
        const idx_t outRowStride = OW * OC;

        const idx_t w = W - 2;  // fixed input col base for this right kernel
        const idx_t out_w =
            w + 1;  // fixed output col we write (shifted right by 1)

        for (idx_t oc = 0; oc < OC; ++oc) {
            const int32_t bias_val = biasBase[oc];

            // first output is at (h=1, out_w, oc) because store uses (h+1)
            int32_t *outPtr = outBase + (1 * outRowStride) + (out_w * OC) + oc;

            for (idx_t h = 0; h < OH - 2; ++h) {
                int32_t partial_sum = bias_val;

                idx_t kh = 0, kw = 0;

                for (idx_t f = 0; f < 6; ++f) {
                    const idx_t in_h = h + kh;
                    const idx_t in_w = w + kw;

                    const idx_t inHWOffset = (in_h * W + in_w) * IC;
                    const idx_t filtHWOffset =
                        filtOcOffset + ((kh * 3) + kw) * IC;

                    vec4xi8_t *inPos = (vec4xi8_t *)(inBase + inHWOffset);
                    vec4xi8_t *kPos = (vec4xi8_t *)(filterBase + filtHWOffset);

                    LOAD_VEC8(x, inPos);  // x0..x7
                    LOAD_VEC8(k, kPos);   // k0..k7
                    SDOT_VEC8(s, x, k);   // s0..s7

                    partial_sum += add_tree_8(s0, s1, s2, s3, s4, s5, s6, s7);

                    kw = (kw == 1) ? 0 : (kw + 1);
                    kh = (kw == 0) ? (kh + 1) : kh;
                }

                *outPtr = partial_sum;   // store at (h+1, out_w, oc)
                outPtr += outRowStride;  // bump to next h (next output row)
            }

            filtOcOffset += filtAdvance;
        }
    }
}

__effcc_rip_exact void conv2d_stride_2_pad_1_impl_32_right(
    // ==== Input (i8, NHWC) ====
    const int8_t *inBase, idx_t W, idx_t H,

    // ==== Filter (i8, [OC,KH,KW,IC]) ====
    const int8_t *filterBase, idx_t OC, idx_t IC,

    // ==== Zero-point ====
    int8_t zp,

    // ==== Bias (i32, length = OC) ====
    const int32_t *biasBase,

    // ==== Output (i8, NHWC) ====
    int32_t *outBase, idx_t OH, idx_t OW) {
    __effcc_ignore_memory_order {
        idx_t filtOcOffset = 0;
        const idx_t filtAdvance = 9 * IC;
        const idx_t outRowStride = OW * OC;

        const idx_t ow = OW - 1;  // fixed output col
        const idx_t w = W - 2;    // fixed input base col for this edge kernel

        for (idx_t oc = 0; oc < OC; ++oc) {
            const int32_t bias_val = biasBase[oc];

            // first output is at (oh=0, ow, oc)
            int32_t *outPtr = outBase + (ow * OC) + oc;

            for (idx_t oh = 0; oh < OH - 1; ++oh) {
                int32_t partial_sum = bias_val;

                const idx_t h = oh << 1;

                idx_t kh = 0, kw = 0;

                for (idx_t f = 0; f < 6; ++f) {
                    const idx_t in_h = h + kh;
                    const idx_t in_w = w + kw;

                    const idx_t inHWOffset = (in_h * W + in_w) * IC;
                    const idx_t filtHWOffset =
                        filtOcOffset + ((kh * 3) + kw) * IC;

                    vec4xi8_t *inPos = (vec4xi8_t *)(inBase + inHWOffset);
                    vec4xi8_t *kPos = (vec4xi8_t *)(filterBase + filtHWOffset);

                    LOAD_VEC8(x, inPos);  // x0..x7
                    LOAD_VEC8(k, kPos);   // k0..k7
                    SDOT_VEC8(s, x, k);   // s0..s7

                    partial_sum += add_tree_8(s0, s1, s2, s3, s4, s5, s6, s7);

                    kw = (kw == 1) ? 0 : (kw + 1);
                    kh = (kw == 0) ? (kh + 1) : kh;
                }

                *outPtr = partial_sum;   // store at (oh, ow, oc)
                outPtr += outRowStride;  // bump to next oh (next output row)
            }

            filtOcOffset += filtAdvance;
        }
    }
}

// top-left-corner
__effcc_rip_exact void conv2d_stride_1_pad_1_impl_32_tl(
    // ==== Input (i8, NHWC) ====
    const int8_t *inBase, idx_t W, idx_t H,

    // ==== Filter (i8, [OC,KH,KW,IC]) ====
    const int8_t *filterBase, idx_t OC, idx_t IC,

    // ==== Zero-point ====
    int8_t zp,

    // ==== Bias (i32, length = OC) ====
    const int32_t *biasBase,

    // ==== Output (i8, NHWC) ====
    int32_t *outBase, idx_t OH, idx_t OW) {
    __effcc_ignore_memory_order {
        idx_t filtOcOffset = 0;
        idx_t filtAdvance = 9 * IC;

        int32_t *outPtr = outBase;  // points to (0,0,oc) as oc increments

        for (idx_t oc = 0; oc < OC; ++oc) {
            const int32_t bias_val = biasBase[oc];
            int32_t partial_sum = bias_val;

            const idx_t w = 0;
            const idx_t h = 0;

            idx_t kh = 0, kw = 0;

            for (idx_t f = 0; f < 4; ++f) {
                const idx_t in_h = h + kh;
                const idx_t in_w = w + kw;

                const idx_t inHWOffset = (in_h * W + in_w) << 5;  // *32

                // taps are shifted (+1,+1) in filter for TL corner
                const idx_t filtHWOffset =
                    filtOcOffset + ((((kh + 1) * 3) + (kw + 1)) << 5);

                vec4xi8_t *inPos = (vec4xi8_t *)(inBase + inHWOffset);
                vec4xi8_t *kPos = (vec4xi8_t *)(filterBase + filtHWOffset);

                LOAD_VEC8(x, inPos);
                LOAD_VEC8(k, kPos);
                SDOT_VEC8(s, x, k);

                partial_sum += add_tree_8(s0, s1, s2, s3, s4, s5, s6, s7);

                kw = (kw == 1) ? 0 : (kw + 1);
                kh = (kw == 0) ? (kh + 1) : kh;
            }

            *outPtr = partial_sum;  // store at (0,0,oc)
            ++outPtr;               // bump to next oc

            filtOcOffset += filtAdvance;
        }
    }
}

// top-right corner
__effcc_rip_exact void conv2d_stride_1_pad_1_impl_32_tr(
    // ==== Input (i8, NHWC) ====
    const int8_t *inBase, idx_t W, idx_t H,

    // ==== Filter (i8, [OC,KH,KW,IC]) ====
    const int8_t *filterBase, idx_t OC, idx_t IC,

    // ==== Zero-point ====
    int8_t zp,

    // ==== Bias (i32, length = OC) ====
    const int32_t *biasBase,

    // ==== Output (i8, NHWC) ====
    int32_t *outBase, idx_t OH, idx_t OW) {
    __effcc_ignore_memory_order {
        idx_t filtOcOffset = 0;
        const idx_t filtAdvance = 9 * IC;

        const idx_t w = W - 2;      // fixed input base col
        const idx_t out_w = w + 1;  // fixed output col (rightmost)

        // output pointer for oc=0 at (h=0, out_w, oc)
        int32_t *outPtr = outBase + (out_w * OC);

        for (idx_t oc = 0; oc < OC; ++oc) {
            const int32_t bias_val = biasBase[oc];
            int32_t partial_sum = bias_val;

            const idx_t h = 0;

            idx_t kh = 0, kw = 0;

            for (idx_t f = 0; f < 4; ++f) {
                const idx_t in_h = h + kh;
                const idx_t in_w = w + kw;

                const idx_t inHWOffset = (in_h * W + in_w) * IC;
                const idx_t filtHWOffset =
                    filtOcOffset + (((kh + 1) * 3) + kw) * IC;

                vec4xi8_t *inPos = (vec4xi8_t *)(inBase + inHWOffset);
                vec4xi8_t *kPos = (vec4xi8_t *)(filterBase + filtHWOffset);

                LOAD_VEC8(x, inPos);
                LOAD_VEC8(k, kPos);
                SDOT_VEC8(s, x, k);

                partial_sum += add_tree_8(s0, s1, s2, s3, s4, s5, s6, s7);

                kw = (kw == 1) ? 0 : (kw + 1);
                kh = (kw == 0) ? (kh + 1) : kh;
            }

            *outPtr = partial_sum;  // store at (0, out_w, oc)
            ++outPtr;               // bump to next channel

            filtOcOffset += filtAdvance;
        }
    }
}

// bottom-left-corner
__effcc_rip_exact void conv2d_stride_1_pad_1_impl_32_bl(
    // ==== Input (i8, NHWC) ====
    const int8_t *inBase, idx_t W, idx_t H,

    // ==== Filter (i8, [OC,KH,KW,IC]) ====
    const int8_t *filterBase, idx_t OC, idx_t IC,

    // ==== Zero-point ====
    int8_t zp,

    // ==== Bias (i32, length = OC) ====
    const int32_t *biasBase,

    // ==== Output (i8, NHWC) ====
    int32_t *outBase, idx_t OH, idx_t OW) {
    __effcc_ignore_memory_order {
        idx_t filtOcOffset = 0;
        const idx_t filtAdvance = 9 * IC;

        const idx_t w = 0;
        const idx_t h = H - 2;

        const idx_t out_h = h + 1;  // bottom output row
        const idx_t outRowStride = OW * OC;

        // output pointer for oc=0 at (out_h, out_w=0, oc)
        int32_t *outPtr = outBase + (out_h * outRowStride);

        for (idx_t oc = 0; oc < OC; ++oc) {
            const int32_t bias_val = biasBase[oc];
            int32_t partial_sum = bias_val;

            idx_t kh = 0, kw = 0;

            for (idx_t f = 0; f < 4; ++f) {
                const idx_t in_h = h + kh;
                const idx_t in_w = w + kw;

                const idx_t inHWOffset = (in_h * W + in_w) * IC;
                const idx_t filtHWOffset =
                    filtOcOffset + ((kh * 3) + (kw + 1)) * IC;

                vec4xi8_t *inPos = (vec4xi8_t *)(inBase + inHWOffset);
                vec4xi8_t *kPos = (vec4xi8_t *)(filterBase + filtHWOffset);

                LOAD_VEC8(x, inPos);
                LOAD_VEC8(k, kPos);
                SDOT_VEC8(s, x, k);

                partial_sum += add_tree_8(s0, s1, s2, s3, s4, s5, s6, s7);

                kw = (kw == 1) ? 0 : (kw + 1);
                kh = (kw == 0) ? (kh + 1) : kh;
            }

            *outPtr = partial_sum;  // store at (out_h, 0, oc)
            ++outPtr;               // bump to next channel

            filtOcOffset += filtAdvance;
        }
    }
}

// bottom-right corner
__effcc_rip_exact void conv2d_stride_1_pad_1_impl_32_br(
    // ==== Input (i8, NHWC) ====
    const int8_t *inBase, idx_t W, idx_t H,

    // ==== Filter (i8, [OC,KH,KW,IC]) ====
    const int8_t *filterBase, idx_t OC, idx_t IC,

    // ==== Zero-point ====
    int8_t zp,

    // ==== Bias (i32, length = OC) ====
    const int32_t *biasBase,

    // ==== Output (i8, NHWC) ====
    int32_t *outBase, idx_t OH, idx_t OW) {
    __effcc_ignore_memory_order {
        idx_t filtOcOffset = 0;
        const idx_t filtAdvance = 9 * IC;

        const idx_t w = W - 2;
        const idx_t h = H - 2;

        const idx_t out_h = h + 1;
        const idx_t out_w = w + 1;

        const idx_t outRowStride = OW * OC;

        // output pointer for oc=0 at (out_h, out_w, oc)
        int32_t *outPtr = outBase + (out_h * outRowStride) + (out_w * OC);

        for (idx_t oc = 0; oc < OC; ++oc) {
            const int32_t bias_val = biasBase[oc];
            int32_t partial_sum = bias_val;

            idx_t kh = 0, kw = 0;

            for (idx_t f = 0; f < 4; ++f) {
                const idx_t in_h = h + kh;
                const idx_t in_w = w + kw;

                const idx_t inHWOffset = (in_h * W + in_w) * IC;
                const idx_t filtHWOffset = filtOcOffset + ((kh * 3) + kw) * IC;

                vec4xi8_t *inPos = (vec4xi8_t *)(inBase + inHWOffset);
                vec4xi8_t *kPos = (vec4xi8_t *)(filterBase + filtHWOffset);

                LOAD_VEC8(x, inPos);
                LOAD_VEC8(k, kPos);
                SDOT_VEC8(s, x, k);

                partial_sum += add_tree_8(s0, s1, s2, s3, s4, s5, s6, s7);

                kw = (kw == 1) ? 0 : (kw + 1);
                kh = (kw == 0) ? (kh + 1) : kh;
            }

            *outPtr = partial_sum;  // store at (out_h, out_w, oc)
            ++outPtr;               // bump to next channel

            filtOcOffset += filtAdvance;
        }
    }
}

__effcc_rip_exact void conv2d_stride_2_pad_1_impl_32_br(
    const int8_t *inBase, idx_t W, idx_t H, const int8_t *filterBase, idx_t OC,
    idx_t IC, int8_t zp, const int32_t *biasBase, int32_t *outBase, idx_t OH,
    idx_t OW) {
    __effcc_ignore_memory_order {
        idx_t filtOcOffset = 0;
        const idx_t filtAdvance = 9 * IC;

        // bottom-right output coordinate
        const idx_t oh = OH - 1;
        const idx_t ow = OW - 1;

        // bottom-right input base (2x2 valid region)
        const idx_t h = H - 2;
        const idx_t w = W - 2;

        const idx_t outRowStride = OW * OC;

        // output pointer for oc=0 at (oh, ow, oc)
        int32_t *outPtr = outBase + (oh * outRowStride) + (ow * OC);

        for (idx_t oc = 0; oc < OC; ++oc) {
            int32_t partial_sum = biasBase[oc];

            idx_t kh = 0, kw = 0;

            // Only 2 rows x 2 cols are valid => 4 taps
            for (idx_t f = 0; f < 4; ++f) {
                const idx_t in_h = h + kh;
                const idx_t in_w = w + kw;

                const idx_t inHWOffset = (in_h * W + in_w) * IC;
                const idx_t filtHWOffset = filtOcOffset + ((kh * 3) + kw) * IC;

                vec4xi8_t *inPos = (vec4xi8_t *)(inBase + inHWOffset);
                vec4xi8_t *kPos = (vec4xi8_t *)(filterBase + filtHWOffset);

                LOAD_VEC8(x, inPos);
                LOAD_VEC8(k, kPos);
                SDOT_VEC8(s, x, k);

                partial_sum += add_tree_8(s0, s1, s2, s3, s4, s5, s6, s7);

                kw = (kw == 1) ? 0 : (kw + 1);
                kh = (kw == 0) ? (kh + 1) : kh;
            }

            *outPtr = partial_sum;  // store at (oh, ow, oc)
            ++outPtr;               // bump to next channel

            filtOcOffset += filtAdvance;
        }
    }
}

__effcc_rip_exact void conv2d_bias_zp_with_hw_indices(
    const int8_t *filter,     // [OC, 3, 3, IC] flattened per OC
    idx_t IC,                 // input channels
    const idx_t *hw_indices,  // list of HW indices (each in [0,8])
    idx_t hw_count,           // number of HW indices
    const int32_t *bias, idx_t bias_dim_oc, int32_t zp, int32_t *bias_w_zp) {
    __effcc_ignore_memory_order {
        const idx_t filter_stride_oc = 9 * IC;

        for (idx_t oc = 0; oc < bias_dim_oc; ++oc) {
            int32_t sumW = 0;

            // sum weights over selected HW taps
            for (idx_t t = 0; t < hw_count; ++t) {
                const idx_t base = hw_indices[t] * IC;
                __effcc_parallel(4) for (idx_t c = 0; c < IC; c += 8) {
                    vec4xi8_t x0 = *(vec4xi8_t *)(filter + base + c);
                    vec4xi8_t x1 = *(vec4xi8_t *)(filter + base + c + 4);

                    int32_t s0 = (int32_t)x0[0];
                    int32_t s1 = (int32_t)x0[1];
                    int32_t s2 = (int32_t)x0[2];
                    int32_t s3 = (int32_t)x0[3];

                    int32_t s4 = (int32_t)x1[0];
                    int32_t s5 = (int32_t)x1[1];
                    int32_t s6 = (int32_t)x1[2];
                    int32_t s7 = (int32_t)x1[3];

                    sumW += add_tree_8(s0, s1, s2, s3, s4, s5, s6, s7);
                }
            }

            // Fold zp into bias: sum((x - zp)*w) + bias
            bias_w_zp[oc] = bias[oc] - zp * sumW;

            filter += filter_stride_oc;  // next OC
        }
    }
}

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

void conv2d_stride_1_pad_1_impl_32(
    // ==== Input (i8, NHWC) ====
    const int8_t *inBase, idx_t W, idx_t H,

    // ==== Filter (i8, [OC,KH,KW,IC]) ====
    const int8_t *filterBase, idx_t OC, idx_t IC, idx_t filter_stride_oc,

    // ==== Zero-point ====
    int8_t zp,

    // ==== Bias (i32, length = OC) ====
    const int32_t *biasBase,

    // ==== Output (i8, NHWC) ====
    int32_t *outBase, idx_t OH, idx_t OW) {
    int32_t bias_w_zp[OC];
    conv2d_bias_zp(filterBase, filter_stride_oc, biasBase, OC, -zp, bias_w_zp);

    conv2d_stride_1_pad_1_impl_32_core(inBase, W, H, filterBase, OC, IC, zp,
                                       bias_w_zp, outBase, OH, OW);

    static const idx_t HW_TOP[6] = {3, 4, 5, 6, 7, 8};
    conv2d_bias_zp_with_hw_indices(filterBase, IC, HW_TOP, 6, biasBase, OC, zp,
                                   bias_w_zp);

    conv2d_stride_1_pad_1_impl_32_top(inBase, W, H, filterBase, OC, IC, zp,
                                      bias_w_zp, outBase, OH, OW);

    static const idx_t HW_BOTTOM[6] = {0, 1, 2, 3, 4, 5};
    conv2d_bias_zp_with_hw_indices(filterBase, IC, HW_BOTTOM, 6, biasBase, OC,
                                   zp, bias_w_zp);

    conv2d_stride_1_pad_1_impl_32_bottom(inBase, W, H, filterBase, OC, IC, zp,
                                         bias_w_zp, outBase, OH, OW);

    static const idx_t HW_LEFT[6] = {1, 2, 4, 5, 7, 8};
    conv2d_bias_zp_with_hw_indices(filterBase, IC, HW_LEFT, 6, biasBase, OC, zp,
                                   bias_w_zp);

    conv2d_stride_1_pad_1_impl_32_left(inBase, W, H, filterBase, OC, IC, zp,
                                       bias_w_zp, outBase, OH, OW);

    static const idx_t HW_RIGHT[6] = {0, 1, 3, 4, 6, 7};
    conv2d_bias_zp_with_hw_indices(filterBase, IC, HW_RIGHT, 6, biasBase, OC,
                                   zp, bias_w_zp);

    conv2d_stride_1_pad_1_impl_32_right(inBase, W, H, filterBase, OC, IC, zp,
                                        bias_w_zp, outBase, OH, OW);

    static const idx_t HW_TL[4] = {4, 5, 7, 8};
    conv2d_bias_zp_with_hw_indices(filterBase, IC, HW_TL, 4, biasBase, OC, zp,
                                   bias_w_zp);

    conv2d_stride_1_pad_1_impl_32_tl(inBase, W, H, filterBase, OC, IC, zp,
                                     bias_w_zp, outBase, OH, OW);

    static const idx_t HW_TR[4] = {3, 4, 6, 7};

    conv2d_bias_zp_with_hw_indices(filterBase, IC, HW_TR, 4, biasBase, OC, zp,
                                   bias_w_zp);

    conv2d_stride_1_pad_1_impl_32_tr(inBase, W, H, filterBase, OC, IC, zp,
                                     bias_w_zp, outBase, OH, OW);

    static const idx_t HW_BL[4] = {1, 2, 4, 5};
    conv2d_bias_zp_with_hw_indices(filterBase, IC, HW_BL, 4, biasBase, OC, zp,
                                   bias_w_zp);

    conv2d_stride_1_pad_1_impl_32_bl(inBase, W, H, filterBase, OC, IC, zp,
                                     bias_w_zp, outBase, OH, OW);

    static const idx_t HW_BR[4] = {0, 1, 3, 4};
    conv2d_bias_zp_with_hw_indices(filterBase, IC, HW_BR, 4, biasBase, OC, zp,
                                   bias_w_zp);

    conv2d_stride_1_pad_1_impl_32_br(inBase, W, H, filterBase, OC, IC, zp,
                                     bias_w_zp, outBase, OH, OW);
}

void conv2d_stride_2_pad_1_impl_32(
    // ==== Input (i8, NHWC) ====
    const int8_t *inBase, idx_t W, idx_t H,

    // ==== Filter (i8, [OC,KH,KW,IC]) ====
    const int8_t *filterBase, idx_t OC, idx_t IC, idx_t filter_stride_oc,

    // ==== Zero-point ====
    int8_t zp,

    // ==== Bias (i32, length = OC) ====
    const int32_t *biasBase,

    // ==== Output (i8, NHWC) ====
    int32_t *outBase, idx_t OH, idx_t OW) {
    int32_t bias_w_zp[OC];
    conv2d_bias_zp(filterBase, filter_stride_oc, biasBase, OC, -zp, bias_w_zp);

    conv2d_stride_2_pad_1_impl_32_core(inBase, W, H, filterBase, OC, IC, zp,
                                       bias_w_zp, outBase, OH, OW);

    static const idx_t HW_BOTTOM[6] = {0, 1, 2, 3, 4, 5};
    conv2d_bias_zp_with_hw_indices(filterBase, IC, HW_BOTTOM, 6, biasBase, OC,
                                   zp, bias_w_zp);

    conv2d_stride_2_pad_1_impl_32_bottom(inBase, W, H, filterBase, OC, IC, zp,
                                         bias_w_zp, outBase, OH, OW);

    static const idx_t HW_RIGHT[6] = {0, 1, 3, 4, 6, 7};
    conv2d_bias_zp_with_hw_indices(filterBase, IC, HW_RIGHT, 6, biasBase, OC,
                                   zp, bias_w_zp);

    conv2d_stride_2_pad_1_impl_32_right(inBase, W, H, filterBase, OC, IC, zp,
                                        bias_w_zp, outBase, OH, OW);

    static const idx_t HW_BR[4] = {0, 1, 3, 4};
    conv2d_bias_zp_with_hw_indices(filterBase, IC, HW_BR, 4, biasBase, OC, zp,
                                   bias_w_zp);

    conv2d_stride_2_pad_1_impl_32_br(inBase, W, H, filterBase, OC, IC, zp,
                                     bias_w_zp, outBase, OH, OW);
}

__efficient__ void pad_input(const int8_t *input, idx_t IC, idx_t IH, idx_t IW,
                             idx_t pad, int8_t zp, int8_t *output) {
    __effcc_ignore_memory_order {
        const vec4xi8_t *input_ptr = (const vec4xi8_t *)input;
        vec4xi8_t *output_ptr = (vec4xi8_t *)(output);
        const vec4xi8_t zp_vec = {zp, zp, zp, zp};
        for (idx_t h = 0; h < IH + 2 * pad; h++) {
            for (idx_t w = 0; w < IW + 2 * pad; w++) {
                for (idx_t c = 0; c < IC; c += 4) {
                    if (h < pad || h >= IH + pad || w < pad || w >= IW + pad) {
                        *output_ptr = zp_vec;
                    } else {
                        *output_ptr = *input_ptr;
                        input_ptr++;
                    }
                    output_ptr++;
                }
            }
        }
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
    if (filter_dim_ic == 32) {
        conv2d_stride_1_pad_1_impl_32(
            in_base + in_offset, in_dim_w, in_dim_h,
            filter_base + filter_offset, filter_dim_oc, filter_dim_ic,
            filter_stride_oc, pad0_ptr[0], bias_base + bias_offset,
            out_base + out_offset, out_dim_h, out_dim_w);

    } else {
        int32_t bias_w_zp[bias_dim_oc];
        conv2d_bias_zp(filter_base + filter_offset, filter_stride_oc,
                       bias_base + bias_offset, bias_dim_oc,
                       -(int32_t)pad0_ptr[0], bias_w_zp);

        conv2d_stride_1_pad_1_impl(in_base + in_offset,
                                   in_dim_w,  // in_dim_h,
                                   filter_base + filter_offset, filter_dim_oc,
                                   filter_dim_ic, pad0_ptr[0], bias_w_zp,
                                   out_base + out_offset,
                                   // temp,
                                   out_dim_h, out_dim_w);
    }
}

__efficient__ void conv2d_stride_2_pad_corner_impl(
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
    vec4xi8_t zp_vec = {zp, zp, zp, zp};
    const int8_t *input_base =
        input;  // we should not adjust the input pointer since there is no top
                // padding for stride-2
    const idx_t IC_div4_2 = 2 * IC_div4;
    const idx_t IW_IC_div4_2 = 2 * IW_IC_div4;

    __effcc_ignore_memory_order {
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
    if (filter_dim_ic == 32) {
        conv2d_stride_2_pad_1_impl_32(
            in_base + in_offset, in_dim_w, in_dim_h,
            filter_base + filter_offset, filter_dim_oc, filter_dim_ic,
            filter_stride_oc, pad0_ptr[0], bias_base + bias_offset,
            out_base + out_offset, out_dim_h, out_dim_w);

    } else {
        int32_t bias_w_zp[bias_dim_oc];
        conv2d_bias_zp(filter_base + filter_offset, filter_stride_oc,
                       bias_base + bias_offset, bias_dim_oc,
                       -(int32_t)pad0_ptr[0], bias_w_zp);
        conv2d_stride_2_pad_corner_impl(
            in_base + in_offset, in_dim_w, filter_base + filter_offset,
            filter_dim_oc, filter_dim_ic, pad0_ptr[0], bias_w_zp,
            out_base + out_offset, out_dim_h, out_dim_w);
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