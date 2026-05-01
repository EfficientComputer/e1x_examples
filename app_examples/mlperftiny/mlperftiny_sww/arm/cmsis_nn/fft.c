// clang-format off <<AUTOBENCH>> ambiq_apollo4p nxp_lpc55s69 renesas_ra8m1 ambiq_apollo510 renesas_ra8p1 clang-format on

#include "dsp/transform_functions.h"
#include "fft4k.h"

arm_cfft_instance_q15 fftInstance;

void fft_init() { arm_cfft_init_q15(&fftInstance, FFT_SIZE); }

void eff_fft4(eff_fft_cpx *src, eff_fft_cpx *dst) {
    // copying src into dst
    for (int i = 0; i < FFT_SIZE; i++) {
        dst[i] = src[i];
    }
    arm_cfft_q15(&fftInstance, (q15_t *)dst, 0, 1);
}
