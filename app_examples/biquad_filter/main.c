#include <eff.h>
#include <eff/profile.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#define DEFINE_SAMPLE_DATA
#include <biquad_filter_input.h.inc>

#define INPUT_SIZE 4096

// Low pass filter centered at 300Hz Q=0.72
float poles_coeffs[2] = {-1.8400889269370049f, 0.8523343179964432f};
float zeros_coeffs[3] = {0.003061347764859629f, 0.006122695529719258f, 0.003061347764859629f};

void biquad_filter(const int16_t* input, int16_t* output, int length, const float p[2], const float z[3]);

int16_t output[INPUT_SIZE];

int main() {
    START_PROFILE_REGION("kernel");
    for (int iter = 0; iter < 100; iter++)
        biquad_filter(sample_input, output, INPUT_SIZE, poles_coeffs, zeros_coeffs);
    END_PROFILE_REGION();

    // Compute error
    int32_t err = 0;
    for (int i = 100; i < INPUT_SIZE; i++) {
        err += abs(output[i] - expected_output[i]);
    }
    
    if ((err / (INPUT_SIZE - 100)) > 660) {
        printf("[biquad_filter] FAIL - Average error exceeded threshold: %d\n", err / (INPUT_SIZE - 100));
        return 1;
    }

    printf("[biquad_filter] PASS\n");
    return 0;
}