#define USE_SUBWORD_SIMD 1

#if USE_SUBWORD_SIMD
typedef char v4sc __attribute__((vector_size(4)));
typedef v4sc data_t;

#define DATA_ZERO {0, 0, 0, 0}
#define DATA_INIT(x) \
    (data_t) { x, x, x, x } // gcc barfs unless you include (data_t)
#else
typedef int data_t;

#define DATA_ZERO 0
#define DATA_INIT(x) x
#endif

// Inner-kernel of dense convolution that assumes a symmetric filter to reduce
// xdata requirements
void dconv5x5sym_inner(const data_t *in_row_zero, data_t *out_row_neg5, int n,
                       int tile_rows_plus_4, data_t f00, data_t f10, data_t f20,
                       data_t f01, data_t f11, data_t f21, data_t f02,
                       data_t f12, data_t f22);
