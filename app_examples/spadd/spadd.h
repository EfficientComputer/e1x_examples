#include <stdint.h>

void spadd(const int32_t *A_vals, const uint32_t *A_rows,
           const uint32_t *A_coords, const int32_t *B_vals,
           const uint32_t *B_rows, const uint32_t *B_coords, int32_t *Z,
           uint32_t rows);
