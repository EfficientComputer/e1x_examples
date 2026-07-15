#include <stdint.h>
typedef int32_t data_t;
#define DATA_ZERO 0
#define DATA_INIT(x) (x)

#define M 64
#define N 64

#define INSTRIDE N
#define OUTSTRIDE N

void dconv3x3(const data_t *in, int m, int n, int mstride, const data_t *filter,
              data_t *out);
