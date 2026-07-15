#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "svec.h"

#ifdef SVEC_DYNAMIC
svec_t *svec(unsigned long nnzs, unsigned long *c, unsigned long *v)
{
    svec_t *r = (svec_t *)malloc(sizeof(svec_t));
    r->nnzs = nnzs;
    r->c = (unsigned long *)malloc(sizeof(unsigned long) * nnzs);
    r->v = (unsigned long *)malloc(sizeof(unsigned long) * nnzs);
    if (c != NULL)
    {
        memcpy(r->c, c, nnzs * sizeof(unsigned long));
        memcpy(r->v, v, nnzs * sizeof(unsigned long));
    }
    return r;
}
#endif

void svec_print(svec_t *a)
{
    printf("[");
    for (int i = 0; i < a->nnzs; i++)
    {
        printf("(%lu %lu)", a->c[i], a->v[i]);
    }
    printf("]\n");
}

svec_t *svec_init_random_vec(svec_t *r)
{
#ifdef SVEC_DYNAMIC
    unsigned long *c =
        (unsigned long *)malloc(MAX_VEC_NNZS * sizeof(unsigned long));
    unsigned long *v =
        (unsigned long *)malloc(MAX_VEC_NNZS * sizeof(unsigned long));
    unsigned long nnz = 0;
    for (int i = 0; i < MAX_COORD; i++)
    {
        if (rand() % NNZ_LIKELIHOOD == 0)
        {
            c[nnz] = i;
            v[nnz] = rand() % CHAR_MAX;
            nnz++;
            if (nnz >= MAX_VEC_NNZS)
            {
                break;
            }
        }
    }

    svec_t *rr = NULL;
    rr = svec(nnz - 1, c, v);
    free(c);
    free(v);
    return rr;

#else

    unsigned long nnz = 0;
    for (int i = 0; i < MAX_COORD; i++)
    {
        if (rand() % NNZ_LIKELIHOOD == 0)
        {
            r->c[nnz] = i;
            r->v[nnz] = rand() % CHAR_MAX;
            nnz++;
            if (nnz >= MAX_VEC_NNZS)
            {
                break;
            }
        }
    }

    r->nnzs = nnz - 1;
    return r;

#endif
}

__efficient__ svec_t *svec_stream_join_mul(svec_t *a, svec_t *b, svec_t *c)
{
    unsigned long aptr = 0;
    unsigned long bptr = 0;
    unsigned long cptr = 0;
    while (aptr < a->nnzs && bptr < b->nnzs)
    {
        unsigned long ac = a->c[aptr];
        unsigned long av = a->v[aptr];
        unsigned long bc = b->c[bptr];
        unsigned long bv = b->v[bptr];

        if (ac == bc)
        {
            c->c[cptr] = ac;
            c->v[cptr] = av * bv;
            aptr++;
            bptr++;
            cptr++;
        }
        else if (ac < bc)
        {
            aptr++;
        }
        else
        {
            bptr++;
        }
    }
    c->nnzs = cptr;
    return c;
}

svec_t *svec_stream_join(svec_t *a, svec_t *b, svec_t *c,
                         unsigned long (*f)(unsigned long, unsigned long))
{
    unsigned long aptr = 0;
    unsigned long bptr = 0;
    unsigned long cptr = 0;
    while (aptr < a->nnzs && bptr < b->nnzs)
    {
        unsigned long ac = a->c[aptr];
        unsigned long av = a->v[aptr];
        unsigned long bc = b->c[bptr];
        unsigned long bv = b->v[bptr];

        if (ac == bc)
        {
            c->c[cptr] = ac;
            c->v[cptr] = f(av, bv);
            aptr++;
            bptr++;
            cptr++;
        }
        else if (ac < bc)
        {
            aptr++;
        }
        else
        {
            bptr++;
        }
    }
    c->nnzs = cptr;
    return c;
}

unsigned long svec_redsum(svec_t *a)
{
    unsigned long r = 0;
    for (int i = 0; i < a->nnzs; i++)
    {
        r += a->v[i];
    }
    return r;
}

unsigned long svec_dot(svec_t *a, svec_t *b, svec_t *c)
{
    svec_stream_join_mul(a, b, c);
    return svec_redsum(c);
}
