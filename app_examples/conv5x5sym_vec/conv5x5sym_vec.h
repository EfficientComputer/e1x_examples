#include <stdint.h>
typedef int32_t data_t;

#define DATA_ZERO 0
#define DATA_INIT(x) (x)

#define M 64
#define N 64

void dconv5x5sym_inner(const data_t *in_row_zero, data_t *out_row_neg5, int n,
                       int tile_rows_plus_4, data_t f00, data_t f10, data_t f20,
                       data_t f01, data_t f11, data_t f21, data_t f02,
                       data_t f12, data_t f22);
