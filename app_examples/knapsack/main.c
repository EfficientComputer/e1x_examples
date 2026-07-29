#include <stdio.h>
#include <stdlib.h>

#include "knapsack.h"

#define NUM_ITERATIONS 1

#define n 5
#define w 400

void init_items(int *weights, int *values, int *count)
{
    weights[0] = 9;
    values[0] = 150;
    count[0] = 1;
    weights[1] = 13;
    values[1] = 35;
    count[1] = 1;
    weights[2] = 153;
    values[2] = 200;
    count[2] = 2;
    weights[3] = 50;
    values[3] = 60;
    count[3] = 2;
    weights[4] = 15;
    values[4] = 60;
    count[4] = 2;
}

int weights[n], values[n], count[n];
int mm[(n + 1) * (w + 1)];
int *m[n + 1];
int s[n];

int main()
{
    int i, tc = 0, tw = 0, tv = 0;
    init_items(weights, values, count);

    for (int i = 0; i < NUM_ITERATIONS; i++)
        knapsack(weights, values, count, w, mm, m, s, n);

    for (i = 0; i < n; i++)
    {
        if (s[i])
        {
            tc += s[i];
            tw += s[i] * weights[i];
            tv += s[i] * values[i];
        }
    }

    if (tc != 6 || tw != 395 || tv != 730)
    {
        printf("[knapsack] FAIL \ncount: %d, weight: %d , value: %d", tc, tw,
               tv);
        return 1;
    }
    else
    {
        printf("[knapsack] PASS\n");
    }

    return 0;
}
