#include <stdbool.h>
#include <stdlib.h>

#define CHECK_MALLOC(x)            \
    if (!(x)) {                    \
        printf("Malloc failed\n"); \
        exit(1);                   \
    }

typedef struct nz_coord {
    int row;
    int col;
    int nz_id;
} nz_coord_t;

typedef struct csr_matrix {
    int *ptrs;
    int *ids;
    int *nz_ids;
} csr_matrix_t;

typedef struct {
    int *mLv2c;
    int *mLc2v;
    int nnz;
    int *F;
    int *B;
    int *codeword_llrs_in;
    int *codeword_llrs_acc;
    int *codeword_decision;
    int *syndrome;
} decode_tmp_allocations_t;

typedef struct {
    csr_matrix_t *row_view;
    csr_matrix_t *col_view;
    nz_coord_t *nz_coords;
    int nnz;
} sparse_matrix_t;

void create_csr_row_view(nz_coord_t *nz_coords, csr_matrix_t *csr, int n_ptrs,
                         int nnz);

void create_csr_col_view(nz_coord_t *nz_coords, csr_matrix_t *csr, int n_ptrs,
                         int nnz);

int create_generator_matrix(int message_length, int codeword_length,
                            int nz_frac_numerator, int nz_frac_denominator,
                            nz_coord_t **generator_nz_coords);

int create_parity_check_matrix(nz_coord_t *generator_nz_coords,
                               int generator_nnz, int message_length,
                               int codeword_length,
                               nz_coord_t **parity_check_nz_coords);
