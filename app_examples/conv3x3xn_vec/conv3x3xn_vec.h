#ifndef CONV3X3XN_VEC_H
#define CONV3X3XN_VEC_H

#include <stdint.h>

#define data_t int32_t
#define DATA_INIT(x) (x)

#define M 16 // number of rows
#define N 16 // number of columns
#define KERNEL_DIM 3
#define NCHANNELS 4

// If this is enabled, we need to include 3 padding columns
// at the beginning of `out` to avoid writing OOB.
#define FASTPATH 1

#if FASTPATH
// INSTRIDE has one padding column to avoid bank conflicts
#define INSTRIDE (N + 1)
// OUTSTRIDE includes the two padding columns at the beginning to
// avoid needing to gate the `out` store
#define OUTSTRIDE (N + KERNEL_DIM - 1)
#else
#define INSTRIDE N
#define OUTSTRIDE N
#endif

#define IS_POWER_OF_TWO(x) ((x) > 0 && ((x) & ((x) - 1)) == 0)

#if !IS_POWER_OF_TWO(N)
#warning "N is not a power of 2, so performance will be poor."
#endif

void dconv_3x3xn(const data_t *in, int m, int n, int instride, int outstride,
                 const data_t *kernels, int nchannels, data_t *out);

#endif
