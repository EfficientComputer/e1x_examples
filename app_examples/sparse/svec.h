#define SVEC_FIXED_VEC
#undef SVEC_DYNAMIC

#define MAX_VEC_NNZS 100
#define SVEC_MAX_NNZS MAX_VEC_NNZS
#define MAX_COORD 100
#define NNZ_LIKELIHOOD 8 /*1 / NNZ_LIKELIHOOD probability that a coord is nz*/

typedef struct svec
{
    unsigned long nnzs;
    unsigned long *c;
    unsigned long *v;
} svec_t;

#ifdef SVEC_DYNAMIC
svec_t *svec();
#endif

svec_t *svec_init_random_vec(svec_t *);
// svec_t *svec_stream_join(svec_t *, svec_t *, svec_t *, unsigned long
// (*)(unsigned long, unsigned long));
svec_t *svec_stream_join_mul(svec_t *, svec_t *, svec_t *);
unsigned long svec_dot(svec_t *, svec_t *, svec_t *);
void svec_print(svec_t *);
