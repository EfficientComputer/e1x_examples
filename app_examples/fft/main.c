#include <eff.h>
#include <eff/profile.h>
#include <stdio.h>

#define DEFINE_SAMPLE_INPUT
#include "fft.h"

fft_cpx out_buf[FFT_SIZE] = {0};

extern int32_t expectedR[FFT_SIZE];
extern int32_t expectedI[FFT_SIZE];

int main() {
    // Run FFT
    START_PROFILE_REGION("kernel");
    for (int iter = 0; iter < 10; iter++)
        fft4((fft_cpx*)sample_input, out_buf);
    END_PROFILE_REGION();

    // Check answer against reference.
    for (int i = 0; FFT_SIZE > i; i++) {
        if (expectedR[i] != out_buf[i].r) {
            printf("[fft4k] FAIL (r:%d) %d\r\n", i, out_buf[i].r);
            return 0;
        }
        if (expectedI[i] != out_buf[i].i) {
            printf("[fft4k] FAIL (i:%d)\r\n", i);
            return 0;
        }
    }

    printf("[fft4k] PASS\r\n");
}