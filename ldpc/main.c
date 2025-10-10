#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>

#include "ldpc.h"

#define MESSAGE_SIZE 336
#define CODEWORD_SIZE 672
#define PARITY_NNZ 2184
#define GENERATOR_NNZ 8736

// This constant is baked into the structure of 802.11 rate 1/2 code. Could vary for other codes.
#define MAX_CHECK_NODE_NNZ 8

// static allocations for 802.11 LDPC test case
extern int H_80211_nnz;
extern int G_80211_nnz;
extern nz_coord_t H_80211[PARITY_NNZ];
extern nz_coord_t G_80211[GENERATOR_NNZ];

extern int generator_col_csr_ptrs_80211[CODEWORD_SIZE+1];
extern int generator_col_csr_ids_80211[GENERATOR_NNZ];
extern int generator_col_csr_nz_ids_80211[GENERATOR_NNZ];

extern int parity_col_csr_ptrs_80211[CODEWORD_SIZE+1];
extern int parity_col_csr_ids_80211[PARITY_NNZ];
extern int parity_col_csr_nz_ids_80211[PARITY_NNZ];

extern int parity_row_csr_ptrs_80211[MESSAGE_SIZE+1];
extern int parity_row_csr_ids_80211[PARITY_NNZ];
extern int parity_row_csr_nz_ids_80211[PARITY_NNZ];

int _test_message[MESSAGE_SIZE];
int _test_codeword[CODEWORD_SIZE];
int _test_codeword_llrs[CODEWORD_SIZE];
int _test_codeword_decision[CODEWORD_SIZE];

// static allocations for LDPC internal workings
int _mLv2c[PARITY_NNZ];
int _mLc2v[PARITY_NNZ];

#ifndef EFF_BLD_HAND_OPTIMIZED
int _F[CODEWORD_SIZE];
int _B[CODEWORD_SIZE];
#else
int _F[CODEWORD_SIZE * MAX_CHECK_NODE_NNZ];
int _B[CODEWORD_SIZE * MAX_CHECK_NODE_NNZ];
#endif

int _codeword_llrs_acc[CODEWORD_SIZE];
int _syndrome[MESSAGE_SIZE];

// Simple test does dynamic allocations
int test() {
    int decoder_max_iterations = 16;
    bool decoder_early_term_possible = true;

    int generator_nnz = G_80211_nnz;

    nz_coord_t *parity_check_nz_coords = H_80211;
    int parity_check_nnz = H_80211_nnz;
    
    csr_matrix_t generator_col_csr;
    generator_col_csr.ptrs = generator_col_csr_ptrs_80211;
    generator_col_csr.ids = generator_col_csr_ids_80211;
    generator_col_csr.nz_ids = generator_col_csr_nz_ids_80211;

    csr_matrix_t parity_row_csr;
    parity_row_csr.ptrs = parity_row_csr_ptrs_80211;
    parity_row_csr.ids = parity_row_csr_ids_80211;
    parity_row_csr.nz_ids = parity_row_csr_nz_ids_80211;
   
    csr_matrix_t parity_col_csr;
    parity_col_csr.ptrs = parity_col_csr_ptrs_80211;
    parity_col_csr.ids = parity_col_csr_ids_80211;
    parity_col_csr.nz_ids = parity_col_csr_nz_ids_80211;
    
    // initialize test message, set some bits
    _test_message[0] = 1;
    _test_message[2] = 1;
    _test_message[9] = 1;
    _test_message[10] = 1;

    for(int i = 0; i < 70; ++i) {
        _test_message[10 + i * 2] = 1;
    }
    
    // encode test_message
    encode(
        &generator_col_csr, 
        _test_message,
        _test_codeword,
        CODEWORD_SIZE);
    
    // initialize test message llrs
    for(int i = 0; i < CODEWORD_SIZE; ++i) {
        // just set the LLRs to -1 or 1 to start
        _test_codeword_llrs[i] = _test_codeword[i] ? -1000 : 1000;
    }

    // introduce a few bit erasures/flips, see if we can error correct
    // bit erasures
    _test_codeword_llrs[0] = 0;
    _test_codeword_llrs[2] = 0;
    _test_codeword_llrs[15] = 0;

    // bit flips
    _test_codeword_llrs[11] = -_test_codeword_llrs[11];
    for(int i = 0; i < 10; ++i) {
        _test_codeword_llrs[15 + i * 40] = -_test_codeword_llrs[15 + i * 40];
    }
    
    sparse_matrix_t parity_matrix;
    parity_matrix.row_view = &parity_row_csr;
    parity_matrix.col_view = &parity_col_csr;
    parity_matrix.nz_coords = parity_check_nz_coords;
    parity_matrix.nnz = parity_check_nnz;
    
    decode_tmp_allocations_t tmp;
    tmp.mLv2c = _mLv2c;
    tmp.mLc2v = _mLc2v;
    tmp.nnz = parity_check_nnz;
    tmp.F = _F;
    tmp.B = _B;
    tmp.codeword_llrs_in = _test_codeword_llrs;
    tmp.codeword_llrs_acc = _codeword_llrs_acc;
    tmp.codeword_decision = _test_codeword_decision;
    tmp.syndrome = _syndrome;

    // decode the test llrs
    int num_iters = decode(
        &parity_matrix,
        &tmp,
        _test_codeword_llrs,
        _test_codeword_decision,
        MESSAGE_SIZE,
        CODEWORD_SIZE,
        MAX_CHECK_NODE_NNZ,
        decoder_max_iterations,
        decoder_early_term_possible);

    for(int i = 0; i < CODEWORD_SIZE; ++i) {
        if(_test_codeword_decision[i] != _test_codeword[i]) {
            printf("[ldpc] FAIL\n");
            return 1;
        }
    }
    printf("LDPC error corrected in %d iterations\n", num_iters);
    printf("[ldpc] PASS\n");

    return 0;
}


int main() {
    test();
}
