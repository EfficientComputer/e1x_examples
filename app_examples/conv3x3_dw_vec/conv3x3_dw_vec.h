#ifndef CONV3X3_DW_VEC_H
#define CONV3X3_DW_VEC_H

#include <stdint.h>

#define data_t int32_t

typedef int8_t data_t_in;   // adjust if your `data_t` is different
typedef int32_t data_t_acc; // accumulator type (matches TOSA i32)
typedef int32_t data_t_out; // output storage type

// Fixed layer spec (can be adjusted):
#define M 48 // number of rows ( height)
#define N 48 // number of columns ( width)
#define NCHANNELS 8
#define OUT_MUL 1
#define KERNEL_DIM 3
#define STRIDE_H 1
#define STRIDE_W 1
#define PAD_TOP 1
#define PAD_LEFT 1
#define PAD_RIGHT 1
#define PAD_BOTTOM 1

#define HEIGHT M
#define WIDTH N

// Derived output sizes
#define OH M //(((M + PAD_TOP + PAD_BOTTOM - KERNEL_DIM) / STRIDE_H) + 1) // 48
#define OW \
    N // (((N + PAD_LEFT + PAD_RIGHT  - KERNEL_DIM) / STRIDE_W) + 1) // 48

void dconv_3x3_dw(const int8_t *in, const int8_t *kernels, int H, int W, int C,
                  int32_t *out);
#endif
