#define DEFINE_SAMPLE_INPUT
#include "fft4k.h"
#include <stdio.h>

fft_cpx out_buf[FFT_SIZE] = {0};

extern int32_t expectedR[FFT_SIZE];
extern int32_t expectedI[FFT_SIZE];

int main() {
    // Run FFT
    fft4((fft_cpx*)sample_input, out_buf);

    // Check answer against reference.
    for (int i = 0; FFT_SIZE > i; i++) {
        if (expectedR[i] != out_buf[i].r) {
            printf("[fft4k] FAIL (r:%d) %d\n", i, out_buf[i].r);
            return 0;
        }
        if (expectedI[i] != out_buf[i].i) {
            printf("[fft4k] FAIL (i:%d)\n", i);
            return 0;
        }
    }

    printf("[fft4k] PASS\n");
}
