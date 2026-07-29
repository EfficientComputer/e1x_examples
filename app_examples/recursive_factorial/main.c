#include <stdio.h>
#include <stdlib.h>

#define NUM_ITERATIONS 1

#define EXPECTED_OUTPUT 3628800

void factorial(int n, int *y);

int main()
{
    int n = 10;
    int fact = 1;

    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
        factorial(n, &fact);

    if (fact != EXPECTED_OUTPUT)
    {
        printf("[recursive_factorial] FAIL -- %d != expected %d\n", fact,
               EXPECTED_OUTPUT);
        return 1;
    }

    printf("[recursive_factorial] PASS\n");
    return 0;
}
