typedef char v4sc __attribute__((vector_size(4)));
typedef v4sc data_t;

#define DATA_ZERO {0, 0, 0, 0}
#define DATA_INIT(x) \
    (data_t) { x, x, x, x }

#define M 64 // number of rows
#define N 64 // number of columns

#define INSTRIDE N
#define OUTSTRIDE N

void dconv3x3(const data_t *in, int m, int n, int mstride, const data_t *filter,
              data_t *out);
