#define FIXED_POINT_BITS 10
#define FIXED_POINT_MUL_RESCALE(x) ((x) >> FIXED_POINT_BITS)
#define FIXED_POINT_DIVIDEND(x) ((int32_t)(x) << FIXED_POINT_BITS)
#define FIXED_POINT_CAST_TO(x) (int32_t)((x) * (1 << FIXED_POINT_BITS))
#define FIXED_POINT_CAST_FROM(x) (int32_t)(x >> FIXED_POINT_BITS)

void cholesky(int *orig, int n, int *chol, int ofs);
