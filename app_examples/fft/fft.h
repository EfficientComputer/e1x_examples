#include <stdint.h>

#define fft_scalar int16_t

typedef struct {
    fft_scalar r;
    fft_scalar i;
} fft_cpx __attribute__((aligned(4)));

void fft(fft_cpx* src, fft_cpx* dst, uint16_t dst_size);
