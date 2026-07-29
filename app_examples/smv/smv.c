#include "smv.h"

#include <stdint.h>
#include <stdlib.h>

#ifdef EFF_BLD_HAND_OPTIMIZED
__attribute__((always_inline)) void smv_loop(const uint32_t *a_indptr,
                                             const uint32_t *a_indices,
                                             const int32_t *a_data,
                                             uint32_t a_nrows, const int32_t *b,
                                             int32_t *z, uint32_t startRow,
                                             uint32_t numLoops)
{
    for (uint32_t row = startRow; row < a_nrows; row += numLoops)
    {
        int w = 0;
        uint32_t endIdx = a_indptr[row + 1];
        for (uint32_t idx = a_indptr[row]; idx != endIdx; idx++)
        {
            uint32_t col = a_indices[idx];
            int a = a_data[idx];
            w += a * b[col];
        }
        z[row] = w;
    }
}

__efficient__ void smv(const uint32_t *a_indptr, const uint32_t *a_indices,
                       const int32_t *a_data, uint32_t a_nrows,
                       const int32_t *b, int32_t *z)
{
    smv_loop(a_indptr, a_indices, a_data, a_nrows, b, z, 0, 6);
    smv_loop(a_indptr, a_indices, a_data, a_nrows, b, z, 1, 6);
    smv_loop(a_indptr, a_indices, a_data, a_nrows, b, z, 2, 6);
    smv_loop(a_indptr, a_indices, a_data, a_nrows, b, z, 3, 6);
    smv_loop(a_indptr, a_indices, a_data, a_nrows, b, z, 4, 6);
    smv_loop(a_indptr, a_indices, a_data, a_nrows, b, z, 5, 6);
}
#else

// matrix A is in sparse CSR format; see scipy.sparse.csr_matrix

__efficient__ void smv(const uint32_t *a_indptr, const uint32_t *a_indices,
                       const int32_t *a_data, uint32_t a_nrows,
                       const int32_t *b, int32_t *z) {
    for (uint32_t row = 0; row < a_nrows; row++) {
        int w = 0;
        for (uint32_t idx = a_indptr[row]; idx != a_indptr[row + 1]; idx++) {
            uint32_t col = a_indices[idx];
            int a = a_data[idx];
            w += a * b[col];
        }
        z[row] = w;
    }
}

#endif
