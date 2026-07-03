#include <stdint.h>

void biquad_filter(const int16_t *input, int16_t *output, int length, const float p[2], const float z[3]);
