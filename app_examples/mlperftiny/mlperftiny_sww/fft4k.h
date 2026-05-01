// clang-format off <<AUTOBENCH>> efficient_e1 efficient_e0 efficient_e1x ambiq_apollo4p nxp_lpc55s69 renesas_ra8m1 ambiq_apollo510 renesas_ra8p1 clang-format on

#pragma once

#ifdef __cplusplus
extern "C" {
#endif
#include "fft4_1024.h"

#define eff_fft_scalar int16_t

typedef struct {
    eff_fft_scalar r;
    eff_fft_scalar i;
} eff_fft_cpx __attribute__((aligned(4)));

void eff_fft4_warmup();

// Only used for competitor ARM devices
void fft_init();

void eff_fft4(eff_fft_cpx* src, eff_fft_cpx* dst);

#ifdef __cplusplus
}
#endif
