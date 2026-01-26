#include <eff.h>
#include <stdio.h>
#include <stdint.h>

typedef int data_t;

#define DATA_ZERO 0
#define DATA_INIT(x) x

void dconv5x5(const data_t *in, int m, int n, const data_t *filter, data_t *out);
