#include <stdint.h>

#include "matrix.h"
#include "qr_decomp.h"

#define NUM_ITERATIONS 1

int main()
{
    setupQRDecomp();

    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
        QRDecomp();

    checkQRDecomp();

    return 0;
}
