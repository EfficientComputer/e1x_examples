#ifndef LDPC_H
#define LDPC_H

#include "ldpc_utils.h"

void encode(
    csr_matrix_t *generator,
    int *message,
    int *codeword,
    int codeword_length);

int decode(
    sparse_matrix_t *parity_matrix,
    decode_tmp_allocations_t *tmp,
    int *codeword_llrs_in,
    int *codeword_decision,
    int message_length,
    int codeword_length,
    int max_check_node_nnz,
    int max_iterations,
    bool early_term_possible);

#endif