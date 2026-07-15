#include <stdint.h>
#include <stdio.h>

#define NUM_ITERATIONS 1

void isqrt(int32_t x, int32_t *result);

int main()
{
    int32_t res;

    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
        isqrt(16, &res);
    if (res != 4)
    {
        printf("isqt(16) -- FAIL\n");
        return 1;
    }

    isqrt(1764, &res);
    if (res != 42)
    {
        printf("isqrt(1764) -- FAIL\n");
        return 1;
    }

    isqrt(20000, &res);
    if (res != 141)
    {
        printf("isqrt(20000) -- FAIL\n");
        return 1;
    }

    printf("[isqrt] Test PASS\n");
    return 0;
}
