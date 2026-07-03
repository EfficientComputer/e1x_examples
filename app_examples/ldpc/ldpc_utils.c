#include "ldpc_utils.h"
#include <stdio.h>

int row_id_comparator(const void *a, const void *b) {
    return (((nz_coord_t *)a)->row < ((nz_coord_t *)b)->row) ? -1 : 1;
}

int col_id_comparator(const void *a, const void *b) {
    return (((nz_coord_t *)a)->col < ((nz_coord_t *)b)->col) ? -1 : 1;
}

void create_csr_view(nz_coord_t *nz_coords, csr_matrix_t *csr, int n_ptrs,
                     int nnz, bool is_row_csr) {
    int *ptrs = (int *)malloc((n_ptrs + 1) * sizeof(int));
    int *ids = (int *)malloc(nnz * sizeof(int));
    int *nz_ids = (int *)malloc(nnz * sizeof(int));
    nz_coord_t *nz_coords_copy = malloc(nnz * sizeof(nz_coord_t));
    if (!ptrs || !ids || !nz_ids || !nz_coords_copy) {
        printf("Failed to allocate csr format\n");
        return;
    }
    csr->ptrs = ptrs;
    csr->ids = ids;
    csr->nz_ids = nz_ids;

    // copy nz coords so we can get info from sorting without
    // modifying underlying data
    for (int i = 0; i < nnz; ++i) {
        nz_coords_copy[i] = nz_coords[i];
    }

    qsort(nz_coords_copy, nnz, sizeof(nz_coord_t),
          is_row_csr ? row_id_comparator : col_id_comparator);

    int nz_idx = 0;
    int last_offset = 0;

    nz_coord_t *nz = &nz_coords_copy[nz_idx];
    int ptr_idx = is_row_csr ? nz->row : nz->col;

    csr->ptrs[0] = 0;

    for (int i = 0; i < n_ptrs; ++i) {
        int n = 0;
        while (ptr_idx == i) {
            n++;
            csr->ids[nz_idx] = is_row_csr ? nz->col : nz->row;
            csr->nz_ids[nz_idx] = nz->nz_id;
            if (nz_idx >= nnz - 1) {
                break;
            }
            nz = &nz_coords_copy[++nz_idx];
            ptr_idx = is_row_csr ? nz->row : nz->col;
        }
        last_offset += n;
        csr->ptrs[i + 1] = last_offset;
    }

    free(nz_coords_copy);
}

void create_csr_col_view(nz_coord_t *nz_coords, csr_matrix_t *csr, int n_ptrs,
                         int nnz) {
    create_csr_view(nz_coords, csr, n_ptrs, nnz, false);
}

void create_csr_row_view(nz_coord_t *nz_coords, csr_matrix_t *csr, int n_ptrs,
                         int nnz) {
    create_csr_view(nz_coords, csr, n_ptrs, nnz, true);
}

void print_csr_matrix(csr_matrix_t *matrix_view, int n_ptrs, int nnz) {
    for (int i = 0; i < n_ptrs + 1; ++i) {
        printf("ptr[%d] = %d\n", i, matrix_view->ptrs[i]);
    }
    for (int i = 0; i < nnz; ++i) {
        printf("ids[%d] = %d, nz_ids[%d] = %d\n", i, matrix_view->ids[i], i,
               matrix_view->nz_ids[i]);
    }
}

int pseudo_rand(int prev) { return (prev * 1000691) & 0x7FFF; }

// creates a standard form generator matrix
int create_generator_matrix(int message_length, int codeword_length,
                            int nz_frac_numerator, int nz_frac_denominator,
                            nz_coord_t **generator_nz_coords) {
    int rows = message_length;
    int cols = codeword_length;

    // this is just an estimate to start
    int nnz_coords = message_length +
                     (nz_frac_numerator * rows * cols) / nz_frac_denominator;
    *generator_nz_coords =
        (nz_coord_t *)malloc(nnz_coords * sizeof(nz_coord_t));
    if (!generator_nz_coords) {
        printf("Couldn't malloc generator nz coords\n");
        return -1;
    }

    // explicitly include message bits
    int nnz = 0;
    for (int i = 0; i < message_length; ++i) {
        (*generator_nz_coords)[nnz].row = i;
        (*generator_nz_coords)[nnz].col = i;
        (*generator_nz_coords)[nnz].nz_id = nnz;
        nnz++;
    }

    // initialize parity checks
    // randomly initialize is_nz
    int pseudo_rand_num = 17;
    for (int i = 0; i < rows; ++i) {
        for (int j = message_length; j < cols; ++j) {
            pseudo_rand_num = pseudo_rand(pseudo_rand_num);
            bool is_nz = pseudo_rand_num <
                         (nz_frac_numerator * (32768 / nz_frac_denominator));
            if (is_nz) {
                (*generator_nz_coords)[nnz].row = i;
                (*generator_nz_coords)[nnz].col = j;
                (*generator_nz_coords)[nnz].nz_id = nnz;
                nnz++;
            }

            // resize if necessary
            if (nnz >= nnz_coords) {
                nnz_coords *= 2;
                *generator_nz_coords = realloc(*generator_nz_coords,
                                               nnz_coords * sizeof(nz_coord_t));
            }
        }
    }
    return nnz;
}

// takes a standard form generator matrix and converts it
// into a parity check matrix
int create_parity_check_matrix(nz_coord_t *generator_nz_coords,
                               int generator_nnz, int message_length,
                               int codeword_length,
                               nz_coord_t **parity_check_nz_coords) {
    int parity_check_nnz = generator_nnz + codeword_length - 2 * message_length;
    *parity_check_nz_coords =
        (nz_coord_t *)malloc(parity_check_nnz * sizeof(nz_coord_t));
    if (!parity_check_nz_coords) {
        printf("Couldn't malloc parity check nz coords\n");
        return -1;
    }

    // first, add the parity check equations by transposing and shifting
    int nnz = 0;
    for (int i = message_length; i < generator_nnz; ++i) {
        int row = generator_nz_coords[i].row;
        int col = generator_nz_coords[i].col;
        (*parity_check_nz_coords)[nnz].col = row;
        (*parity_check_nz_coords)[nnz].row = col - message_length;
        (*parity_check_nz_coords)[nnz].nz_id = nnz;
        nnz++;
    }

    // then, add a (codeword_length - message_length) x (codeword_length -
    // message_length) identity matrix
    int n_parity_bits = codeword_length - message_length;
    for (int i = 0; i < n_parity_bits; ++i) {
        (*parity_check_nz_coords)[nnz].row = i;
        (*parity_check_nz_coords)[nnz].col = message_length + i;
        (*parity_check_nz_coords)[nnz].nz_id = nnz;
        nnz++;
    }

    return nnz;
}
