#include <stdint.h>

__efficient__ void isqrt(int32_t x, int32_t *result)
{
    int32_t q = 1, r = 0;
    while (q <= x)
    {
        q <<= 2;
    }
    while (q > 1)
    {
        int32_t t;
        q >>= 2;
        t = x - r - q;
        r >>= 1;
        if (t >= 0)
        {
            x = t;
            r += q;
        }
    }
    *result = r;
}
