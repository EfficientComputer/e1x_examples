#ifndef LDPC_UTILS_H
#define LDPC_UTILS_H

#include <stdlib.h>
#include <stdbool.h>

#define CHECK_MALLOC(x)        \
if(!(x)) {                     \
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

#endif