#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "svec.h"

#ifdef SVEC_FIXED_VEC
#define V1LEN 7
unsigned long v1_index[V1LEN] = {19, 20, 35, 41, 60, 76, 79};
unsigned long v1_val[V1LEN] = {74, 83, 118, 0, 45, 40, 11};

#define V2LEN 11
unsigned long v2_index[V2LEN] = {5, 19, 21, 28, 35, 48, 59, 64, 82, 83, 89};
unsigned long v2_val[V2LEN] = {47, 33, 16, 32, 120, 14, 119, 66, 96, 97, 118};

#define RESLEN 2
unsigned res_index[RESLEN] = {19, 35};
unsigned res_val[RESLEN] = {2442, 14160};
#endif

#ifndef SVEC_DYNAMIC
svec_t aa;
unsigned long a_c[SVEC_MAX_NNZS];
unsigned long a_v[SVEC_MAX_NNZS];

svec_t bb;
unsigned long b_c[SVEC_MAX_NNZS];
unsigned long b_v[SVEC_MAX_NNZS];

svec_t cc;
unsigned long c_c[SVEC_MAX_NNZS];
unsigned long c_v[SVEC_MAX_NNZS];
#endif

__efficient__ void f() {}

int main()
{

#ifdef SVEC_DYNAMIC
    svec_t *a = svec_init_random_vec(NULL);
    svec_t *b = svec_init_random_vec(NULL);
    svec_t *c = svec(a->nnzs > b->nnzs ? a->nnzs : b->nnzs, NULL, NULL);
#else
    aa.c = a_c;
    aa.v = a_v;
    bb.c = b_c;
    bb.v = b_v;
    cc.c = c_c;
    cc.v = c_v;
    svec_t *a = &aa;
    svec_t *b = &bb;
    svec_t *c = &cc;

#ifdef SVEC_FIXED_VEC
    for (int i = 0; i < V1LEN; i++)
    {
        a->c[i] = v1_index[i];
        a->v[i] = v1_val[i];
    }
    a->nnzs = V1LEN;

    for (int i = 0; i < V2LEN; i++)
    {
        b->c[i] = v2_index[i];
        b->v[i] = v2_val[i];
    }
    b->nnzs = V2LEN;
#else
    svec_init_random_vec(a);
    svec_init_random_vec(b);
#endif

#endif

    svec_print(a);
    svec_print(b);

    svec_stream_join_mul(a, b, c);
    svec_print(c);

#ifdef SVEC_FIXED_VEC
    for (int i = 0; i < RESLEN; i++)
    {
        if (c->c[i] != res_index[i] || c->v[i] != res_val[i])
        {
            printf("FAIL\n");
            return 0xdead;
        }
    }
    printf("PASS\n");
    return 0;
#endif
}
