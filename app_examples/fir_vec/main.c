#include <stdio.h>

#include "fir_vec.h"

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

    // Verify all N outputs against inline scalar reference (not just o[0]).
    for (int i = 0; i < N; i++)
    {
        int ref = 0;
        for (int j = 0; j < W; j++)
            ref += w[j] * x[i + j];
        if (o[i] != ref)
        {
            printf("[FIR] FAIL at i=%d: got %d expected %d\n", i, o[i], ref);
            return 1;
        }
    }

    printf("[FIR] PASS\n");
    return 0;
}
