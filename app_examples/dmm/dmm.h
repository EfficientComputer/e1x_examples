#include <stdint.h>

void dmm(int32_t *a, // n x m
         int32_t *b, // m x o
         int32_t *z, // n x o
         uint32_t n, uint32_t m, uint32_t o);

void dmm_reference(
    int32_t *a, // n x m
    int32_t *b, // m x o
    int32_t *z, // n x o
    uint32_t n, uint32_t m, uint32_t o);
