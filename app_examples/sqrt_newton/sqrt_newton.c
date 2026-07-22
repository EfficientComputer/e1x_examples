#include <math.h>
#include <stdio.h>

// Until we have softfloat support, we only use float in the fabric,
// and multiply the result by a million to get 6 or so digits of precision
// as an integer.

#ifdef __EFFCC__
__effcc_rip_exact
#endif
    float sqrt_newton(float a)
{
    if (a == 0)
        return 0;

    if (a < 0)
        return -1;

    float x = a;
    float approximation = x / 2;

    for (int i = 0; i < 5; i++)
    {
        float fx = approximation * approximation - x;
        float dx = 2 * approximation;
        approximation = approximation - fx / dx;
    }

    // since n22 can't handle floats, use this to see more
    // digits of precision.
    return approximation;
}

#define SQRT_QTY 10

static float actual_sqrt[] = {
    0,
    1.000000,
    1.414213,
    1.732050,
    2.000000,
    2.236068,
    2.449489,
    2.645751,
    2.828427,
    3.000000,
    3.162277,
};

int float_equality(float a, float b) { return fabsf(a - b) < 0.000001; }

int main()
{
    printf("Calculating square roots\n");

    int results[SQRT_QTY] = {0};

    int pass = 1;

    for (int i = 0; i < SQRT_QTY; i++)
    {
        float in = (float)i;
        float result = sqrt_newton(i);
        results[i] = result;
        printf("Sqrt of %f is %f ", in, result);
        if (float_equality(actual_sqrt[i], result))
        {
            printf("\n");
        }
        else
        {
            printf("but expected %f\n", actual_sqrt[i]);
            pass = 0;
        }
    }

    if (pass)
    {
        printf("[PASS] Sqrt good\n");
    }
    else
    {
        printf("[FAIL] Sqrt bad\n");
        return 1;
    }

    return 0;
}
