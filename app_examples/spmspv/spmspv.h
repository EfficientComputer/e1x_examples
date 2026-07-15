#include <stdint.h>

void spmspv(const uint32_t *a_indptr, const uint32_t *a_indices,
            const int32_t *a_data, uint32_t a_nrows, const uint32_t *x_indices,
            const int32_t *x_data, const int32_t xnnz, int32_t *z);
