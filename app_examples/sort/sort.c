#include "sort.h"

#include <stdbool.h>
#include <stdint.h>

__efficient__ void merge_sort(int32_t *arr, int32_t *temp, int32_t num)
{
    int32_t *dst = temp;
    int32_t *src = arr;

    // k is the level of "recursion" - building up from sorting
    // 2 element arrays until it sorts the entire array.
    for (uint32_t k = 1; k < num; k *= 2)
    {
        // for each array slice
        for (uint32_t left = 0; left < num - k; left += k * 2)
        {
            uint32_t right = left + k;
            uint32_t right_end = ((right + k) > num) ? num : right + k;

            // interleave the left and right sides
            // note that they each are already sorted
            uint32_t m = left, i = left, j = right;
            int32_t a = src[i];
            int32_t b = src[j];
            while (i < right && j < right_end)
            {
                bool decider = a < b;
                dst[m++] = (decider) ? a : b;
                i += (uint32_t)decider;
                j += (uint32_t)!decider;
                int32_t idx = (decider) ? i : j;
                int32_t t = src[idx];
                a = (decider) ? t : a;
                b = (decider) ? b : t;
            }

            uint32_t t = (i < right) ? i : j;
            uint32_t end = (i < right) ? right : right_end;
            for (; t < end; t++)
            {
                dst[m++] = src[t];
            }
        }

        // Swap buffers
        int32_t *tmp = dst;
        dst = src;
        src = tmp;
    }

    uint32_t bound = (src != arr) ? num : 0;
    for (int i = 0; i < bound; i++)
    {
        arr[i] = src[i];
    }
}
