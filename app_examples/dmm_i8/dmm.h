#include <stdint.h>

void dmm(uint32_t *a, // n x m
         uint32_t *b, // m x o
         uint32_t *z, // n x o
         int32_t *zUnscaled, uint32_t n, uint32_t m, uint32_t o,
         int32_t multiplier, int shift);
