#include <stdint.h>
typedef int32_t data_t;

#define DATA_INIT(x) (x)

void dconv5x5(const data_t *in, int m, int n, const data_t *filter,
              data_t *out);
