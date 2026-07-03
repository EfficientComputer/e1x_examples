#include <stdio.h>
#include "fft.h"

// If changed, the CMakeLists.txt file will need to be updated so
// generate_fft.py follows the new size
#define FFT_SIZE 4096
#define abs(x) ((x) < 0 ? -(x) : (x))

// Generated expected data
extern const int32_t expectedR[FFT_SIZE];
extern const int32_t expectedI[FFT_SIZE];
extern const uint32_t sample_input[FFT_SIZE];

fft_cpx out_buf[FFT_SIZE] = {0};

int main() {
    // Run FFT
    fft((fft_cpx*)sample_input, out_buf, FFT_SIZE);

    // Check answer against reference.
    for (int i = 0; FFT_SIZE > i; i++) {
        if (abs(expectedR[i] - out_buf[i].r) > 10 ||
            abs(expectedI[i] - out_buf[i].i) > 10) {
            printf("[fft] FAIL (%d) expected=(%d,%d) got=(%d,%d)\n", i,
                   expectedR[i], expectedI[i], out_buf[i].r, out_buf[i].i);
            return 0;
        }
    }

    printf("[fft] PASS\n");
    return 0;
}
