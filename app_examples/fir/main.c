#include <stdio.h>
#include "fir.h"

#define NUM_ITERATIONS 1

int main()
{
    int w[W], x[N + W - 1], o[N];
    for (int i = 0; i < N + W - 1; i++)
    {
        x[i] = i;
    }
    for (int i = 0; i < W; i++)
    {
        w[i] = i - 8;
    }
    for (int i = 0; i < N; i++)
    {
        o[i] = 0;
    }

    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
        fir(x, w, o);

    printf("[FIR] Expected: 280, Actual: %d -- ", o[0]);
    if (o[0] == 280)
    {
        printf("PASS\n");
    }
    else
    {
        printf("FAIL\n");
        return 1;
    }
    return 0;
}
