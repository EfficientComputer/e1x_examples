#include <eff.h>
#include <stdio.h>
#include <stdint.h>

#define FIXED_ROUND(x) \
    (((x) + (1 << (FRACBITS - 1))) >> FRACBITS)

#define SAMP_MAX 32767
#define FRACBITS 15

__efficient__
void biquad_filter_q15(const int16_t* input, int16_t* output, int length,
                      const int32_t p1, const int32_t p2,
                      const int32_t z0, const int32_t z1, const int32_t z2) {

    int32_t x1 = 0;
    int32_t x2 = 0;
    int32_t y1 = 0;
    int32_t y2 = 0;
  
    for (int n=0; n < length; n +=4)  {

           // ---- sample n ---- 
        int32_t x = input[n + 0];
        int32_t acc = x * z0 + x1 * z1 + x2 * z2
                    - y1 * p1 - y2 * p2;
        int32_t y = FIXED_ROUND(acc);
        output[n + 0] = (int16_t)y;
        // update state
        x2 = x1; x1 = x;
        y2 = y1; y1 = y;

        // ---- sample n+1 ----
        x = input[n + 1];
        acc = x * z0 + x1 * z1 + x2 * z2
                    - y1 * p1 - y2 * p2;
        y = FIXED_ROUND(acc);
        output[n + 1] = (int16_t)y;
        x2 = x1; x1 = x;
        y2 = y1; y1 = y;

        // ---- sample n+2 ----
        x = input[n + 2];
        acc = x * z0 + x1 * z1 + x2 * z2
                    - y1 * p1 - y2 * p2;
        y = FIXED_ROUND(acc);
        output[n + 2] = (int16_t)y;
        x2 = x1; x1 = x;
        y2 = y1; y1 = y;

        // ---- sample n+3 ----
        x = input[n + 3];
        acc = x * z0 + x1 * z1 + x2 * z2
                    - y1 * p1 - y2 * p2;
        y = FIXED_ROUND(acc);
        output[n + 3] = (int16_t)y;
        x2 = x1; x1 = x;
        y2 = y1; y1 = y;  
    } 
}

void biquad_filter(const int16_t* input, int16_t* output, int length, const float p[2], const float z[3]) {
    biquad_filter_q15(input, output, length,
                     (int32_t)(p[0] * SAMP_MAX),
                     (int32_t)(p[1] * SAMP_MAX),
                     (int32_t)(z[0] * SAMP_MAX),
                     (int32_t)(z[1] * SAMP_MAX),
                     (int32_t)(z[2] * SAMP_MAX));
}