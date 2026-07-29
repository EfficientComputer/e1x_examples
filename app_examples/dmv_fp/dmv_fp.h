#include <stdint.h>

#ifdef __EFFCC__
typedef __fp16 float_type;
#else
typedef float float_type;
#endif

#define DATA_INIT(x) ((data_t)(x))

void dmv_fp(float_type *a, float_type *b, float_type *z, uint32_t n, uint32_t m, uint32_t nstride);
