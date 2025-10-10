#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "ldpc.h"

#include <effintrinsic.h>

#ifdef EFF_BLD_HAND_OPTIMIZED
__effcc_rip
#endif
void binary_csr_multiply(
    csr_matrix_t *mat, 
    int *vec_in, 
    int *vec_out,
    int vec_out_len) {
    
    __effcc_ignore_memory_order {
    for(int j = 0; j < vec_out_len; ++j) {
        int offset_start = mat->ptrs[j];
        int offset_end = mat->ptrs[j+1];

        int acc = 0;
        for(int i = offset_start; i < offset_end; ++i) {
            int id = mat->ids[i];
            // add mod 2
            acc ^= vec_in[id];
        }
        vec_out[j] = acc;
    }
    }
}

void encode(
    csr_matrix_t *generator,
    int *message,
    int *codeword,
    int codeword_length) {
    binary_csr_multiply(generator, message, codeword, codeword_length);
}

bool is_codeword(
    csr_matrix_t *parity_check,
    int *codeword_decision,
    int *syndrome,
    int syndrome_length) {
    
    binary_csr_multiply(
        parity_check, 
        codeword_decision, 
        syndrome, 
        syndrome_length);
    
    for(int i = 0; i < syndrome_length; ++i) {
        if(syndrome[i] != 0) {
            return false;
        }
    }
    return true;
}

int sign(int x) {
    return x < 0 ? -1 : 1;
}

int min(int x, int y) {
    return x < y ? x : y;
}

int mCNApprox(int x, int y) {
    return sign(x) * sign(y) * min(abs(x), abs(y));
}

#ifndef EFF_BLD_HAND_OPTIMIZED
__effcc_rip
int decode_(
    sparse_matrix_t *parity_matrix,
    decode_tmp_allocations_t *tmp,
    int message_length,
    int codeword_length,
    int max_check_node_nnz,
    int max_iterations,
    bool early_term_possible) {
    
    csr_matrix_t *parity_row_csr_view = parity_matrix->row_view;
    csr_matrix_t *parity_col_csr_view = parity_matrix->col_view;
    nz_coord_t *parity_nz_coords = parity_matrix->nz_coords;

    int *mLv2c = tmp->mLv2c; // messages from message nodes to check nodes
    int *mLc2v = tmp->mLc2v; // messages from check nodes to message nodes
    int nnz = tmp->nnz;
    int *F = tmp->F; // forward scan of min-sum approximation function
    int *B = tmp->B; // backward scan of min-sum approximation function
    int *codeword_llrs_in = tmp->codeword_llrs_in;
    int *codeword_llrs_acc = tmp->codeword_llrs_acc;
    int *codeword_decision = tmp->codeword_decision;
    int *syndrome = tmp->syndrome;
    
    for(int i = 0; i < nnz; ++i) {
        mLv2c[i] = codeword_llrs_in[parity_nz_coords[i].col];
    }

    int syndrome_length = codeword_length - message_length;
    
    unsigned I = 0;
    while (I < max_iterations) {
        // CN processing
        int n_check_nodes = codeword_length - message_length;
        for (int i = 0; i < n_check_nodes; ++i)
        {
            int check_node_nz_start = parity_row_csr_view->ptrs[i];
            int check_node_nz_end = parity_row_csr_view->ptrs[i+1];
            int check_node_nnz = check_node_nz_end - check_node_nz_start;
            int *nz_ids = parity_row_csr_view->nz_ids + check_node_nz_start;
            
            F[0] = mLv2c[nz_ids[0]];
            B[check_node_nnz - 1] = mLv2c[nz_ids[check_node_nnz - 1]];
            for (int j = 1; j < check_node_nnz; ++j)
            {
                F[j] = mCNApprox(F[j - 1], mLv2c[nz_ids[j]]);
                B[check_node_nnz - 1 - j] = mCNApprox(
                    B[check_node_nnz - j], mLv2c[nz_ids[check_node_nnz - j - 1]]);
            }
            mLc2v[nz_ids[0]] = B[1];
            mLc2v[nz_ids[check_node_nnz - 1]] = F[check_node_nnz - 2];
            for (int j = 1; j < check_node_nnz - 1; ++j)
            {
                mLc2v[nz_ids[j]] = mCNApprox(F[j - 1], B[j + 1]);
            }
        }
        
        // VN processing and app calc
        for (int i = 0; i < codeword_length; ++i) // only transmitted bits
        {
            codeword_llrs_acc[i] = codeword_llrs_in[i];
            
            int msg_node_nz_start = parity_col_csr_view->ptrs[i];
            int msg_node_nz_end = parity_col_csr_view->ptrs[i+1];
            int msg_node_nnz = msg_node_nz_end - msg_node_nz_start;
            int *nz_ids = parity_col_csr_view->nz_ids + msg_node_nz_start;
            for (int j = 0; j < msg_node_nnz; ++j) 
            {
                codeword_llrs_acc[i] += mLc2v[nz_ids[j]];
            }

            codeword_decision[i] = (codeword_llrs_acc[i] <= 0); // approx decision on ith bits

            for (int j = 0; j < msg_node_nnz; ++j)
            {
                mLv2c[nz_ids[j]] = codeword_llrs_acc[i] - mLc2v[nz_ids[j]];
            }
        }

        if (early_term_possible)
        {
            if (is_codeword(
                parity_row_csr_view, 
                codeword_decision,
                syndrome,
                syndrome_length))
            {
                break;
            }
        }

        ++I;
    }

    return I;
}
#else

__effcc_rip 
void initialize_check_node_estimates(
    int * restrict mLv2c,
    int * restrict codeword_llrs_in,
    nz_coord_t * restrict parity_nz_coords,
    int nnz) {
    for(int i = 0; i < nnz; ++i) {
        mLv2c[i] = codeword_llrs_in[parity_nz_coords[i].col];
    }
}

__effcc_rip
void update_message_nodes0(
    int codeword_length,
    int *codeword_llrs_acc,
    int *codeword_llrs_in,
    int *codeword_decision,
    csr_matrix_t *parity_col_csr_view,
    int *mLc2v) {
    __effcc_ignore_memory_order {
    for (int i = 0; i < codeword_length; ++i) // only transmitted bits
    {
        int acc = codeword_llrs_in[i];
        
        int msg_node_nz_start = parity_col_csr_view->ptrs[i];
        int msg_node_nz_end = parity_col_csr_view->ptrs[i+1];
        int msg_node_nnz = msg_node_nz_end - msg_node_nz_start;
        int *nz_ids = parity_col_csr_view->nz_ids + msg_node_nz_start;
        // msg_node_nnz will be small, like O(10)! The LD in LDPC means low-density.
        // Much more efficient to parallelize the outer loop.
        __effcc_parallel(1)
        for (int j = 0; j < msg_node_nnz; ++j) 
        {
            acc += mLc2v[nz_ids[j]];
        }

        codeword_decision[i] = (acc <= 0); // approx decision on ith bits
        codeword_llrs_acc[i] = acc;
    }
    }
}

__effcc_rip
void update_message_nodes1(
    int codeword_length,
    int *codeword_llrs_acc,
    int *codeword_llrs_in,
    int *codeword_decision,
    csr_matrix_t *parity_col_csr_view,
    int *mLc2v,
    int *mLv2c) {
    __effcc_ignore_memory_order {
    for (int i = 0; i < codeword_length; ++i) {
        int msg_node_nz_start = parity_col_csr_view->ptrs[i];
        int msg_node_nz_end = parity_col_csr_view->ptrs[i+1];
        int msg_node_nnz = msg_node_nz_end - msg_node_nz_start;
        int *nz_ids = parity_col_csr_view->nz_ids + msg_node_nz_start;
        __effcc_parallel(1)
        for (int j = 0; j < msg_node_nnz; ++j) {
            mLv2c[nz_ids[j]] = codeword_llrs_acc[i] - mLc2v[nz_ids[j]];
        }
    }
    }
}

__effcc_rip
void compute_f_and_b(
    int n_check_nodes,
    int max_check_node_nnz,
    csr_matrix_t *parity_row_csr_view,
    int *F,
    int *B,
    int *mLv2c) {
    __effcc_ignore_memory_order {
    __effcc_parallel(3)
    for (int i = 0; i < n_check_nodes; ++i){
        int check_node_nz_start = parity_row_csr_view->ptrs[i];
        int check_node_nz_end = parity_row_csr_view->ptrs[i+1];
        int check_node_nnz = check_node_nz_end - check_node_nz_start;
        int *nz_ids = parity_row_csr_view->nz_ids + check_node_nz_start;

        int ftemp = mLv2c[nz_ids[0]];
        int fprevsign = ftemp > 0 ? 1 : -1;
        int fprevmagnitude = fprevsign > 0 ? ftemp : -ftemp;
        
        int btemp = mLv2c[nz_ids[check_node_nnz - 1]];
        int bprevsign = btemp > 0 ? 1 : -1;
        int bprevmagnitude = bprevsign > 0 ? btemp : -btemp;
        
        F[i * max_check_node_nnz + 0] = ftemp;
        B[i * max_check_node_nnz + check_node_nnz - 1] = btemp;
        
        // There's three tricks at play here:
        // 1. Allocate separate F and B accumulators for each check node. This way
        //    we can process all the check nodes in parallel at the expense of some
        //    extra memory usage.
        // 2. Make carry depedencies across inner loop iterations explicit. This allows
        //    us to wrap the entire function with ignore_memory_order, but still preserve
        //    logical ordering in the inner loop.
        // 3. Decompose the mCNApprox min-sum function into its constituent parts. Notice
        //    that even though the recurrence relation for F is 
        //    F[j] = mCNApprox(F[j - 1], mLv2c[nz_ids[j]]), we can actually start computing 
        //    F[j] without fully calculating F[j - 1]. Specifically, the sign of F[j] can 
        //    be computed from that of F[j - 1] and the magnitude of F[j] can be computed 
        //    from that of F[j - 1]. Then, F[j] is simply sign(F[j]) * magnitude(F[j]).
        //    The key observation is both the sign and magnitude recurrences can be implemented 
        //    by carrying F[j - 1] only a single hop back to be compared with F[j]. This reduces
        //    the initiation interval of the inner loops to 1 or 2. They're 4 or 5 without this 
        //    optimization.
        for (int j = 1; j < check_node_nnz; ++j) {
            int mnext = mLv2c[nz_ids[j]];
            int msign = mnext > 0 ? 1 : -1;
            int fnextsign = mnext > 0 ? fprevsign : -fprevsign;
            int mnextabs = msign * mnext;
            int fnextmagnitude = mnextabs < fprevmagnitude ? mnextabs : fprevmagnitude;
            int fnext = fnextsign > 0 ? fnextmagnitude : -fnextmagnitude;

            int mnext2 = mLv2c[nz_ids[check_node_nnz - j - 1]];
            int msign2 = mnext2 > 0 ? 1 : -1;
            int bnextsign = mnext2 > 0 ? bprevsign : -bprevsign;
            int mnextabs2 = msign2 * mnext2;
            int bnextmagnitude = mnextabs2 < bprevmagnitude ? mnextabs2 : bprevmagnitude;
            int bnext = bnextsign > 0 ? bnextmagnitude : -bnextmagnitude;
            F[i * max_check_node_nnz + j] = fnext;
            B[i * max_check_node_nnz + check_node_nnz - 1 - j] = bnext;

            fprevsign = fnextsign;
            fprevmagnitude = fnextmagnitude;
            bprevsign = bnextsign;
            bprevmagnitude = bnextmagnitude;
        }
    }
    }
}

__effcc_rip
void update_c2v(
    int n_check_nodes,
    int max_check_node_nnz,
    csr_matrix_t* parity_row_csr_view,
    int *F,
    int *B,
    int *mLc2v) {
    __effcc_ignore_memory_order {
    for (int i = 0; i < n_check_nodes; ++i) {
        int check_node_nz_start = parity_row_csr_view->ptrs[i];
        int check_node_nz_end = parity_row_csr_view->ptrs[i+1];
        int check_node_nnz = check_node_nz_end - check_node_nz_start;
        int *nz_ids = parity_row_csr_view->nz_ids + check_node_nz_start;
        
        mLc2v[nz_ids[0]] = B[i * max_check_node_nnz + 1];
        mLc2v[nz_ids[check_node_nnz - 1]] = F[i * max_check_node_nnz + check_node_nnz - 2];
        for (int j = 1; j < check_node_nnz - 1; ++j)
        {
            mLc2v[nz_ids[j]] = 
            mCNApprox(F[i * max_check_node_nnz + j - 1], B[i * max_check_node_nnz + j + 1]);
        }
    }
    }
}

int decode_(
    sparse_matrix_t *parity_matrix,
    decode_tmp_allocations_t *tmp,
    int message_length,
    int codeword_length,
    int max_check_node_nnz,
    int max_iterations,
    bool early_term_possible) {
    
    csr_matrix_t *parity_row_csr_view = parity_matrix->row_view;
    csr_matrix_t *parity_col_csr_view = parity_matrix->col_view;
    nz_coord_t *parity_nz_coords = parity_matrix->nz_coords;

    int *mLv2c = tmp->mLv2c; // messages from message nodes to check nodes
    int *mLc2v = tmp->mLc2v; // messages from check nodes to message nodes
    int nnz = tmp->nnz;
    int *F = tmp->F; // forward scan of min-sum approximation function
    int *B = tmp->B; // backward scan of min-sum approximation function
    int *codeword_llrs_in = tmp->codeword_llrs_in;
    int *codeword_llrs_acc = tmp->codeword_llrs_acc;
    int *codeword_decision = tmp->codeword_decision;
    int *syndrome = tmp->syndrome;
    
    initialize_check_node_estimates(
        mLv2c, codeword_llrs_in, parity_nz_coords, nnz);
    
    int syndrome_length = codeword_length - message_length;

    unsigned I = 0;
    while (I < max_iterations) {
        int n_check_nodes = codeword_length - message_length;
        
        compute_f_and_b(
            n_check_nodes,
            max_check_node_nnz,
            parity_row_csr_view,
            F,
            B,
            mLv2c);

        update_c2v(
            n_check_nodes,
            max_check_node_nnz,
            parity_row_csr_view,
            F,
            B,
            mLc2v);
        
        update_message_nodes0(
            codeword_length,
            codeword_llrs_acc,
            codeword_llrs_in,
            codeword_decision,
            parity_col_csr_view,
            mLc2v);
            
        update_message_nodes1(
            codeword_length,
            codeword_llrs_acc,
            codeword_llrs_in,
            codeword_decision,
            parity_col_csr_view,
            mLc2v,
            mLv2c);
        
        if (early_term_possible)
        {
            if (is_codeword(
                parity_row_csr_view, 
                codeword_decision,
                syndrome,
                syndrome_length))
            {
                break;
            }
        }

        ++I;
    }

    return I;
}
#endif

int decode(
    sparse_matrix_t *parity_matrix,
    decode_tmp_allocations_t *tmp,
    int *codeword_llrs_in,
    int *codeword_decision,
    int message_length,
    int codeword_length,
    int max_check_node_nnz,
    int max_iterations,
    bool early_term_possible) {

    int iters_taken;
    iters_taken = decode_(
        parity_matrix,
        tmp,
        message_length,
        codeword_length,
        max_check_node_nnz,
        max_iterations,
        early_term_possible);
    
    return iters_taken;
}
