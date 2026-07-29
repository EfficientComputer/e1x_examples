#ifdef EFF_BLD_HAND_OPTIMIZED
__attribute__((always_inline)) void knapsack_inner_loop(int weight, int value,
                                                        int count, int **m,
                                                        int i, int w,
                                                        int loopIdx,
                                                        int numLoops)
{
    for (int j = loopIdx; j < w; j += numLoops)
    {
        int max = m[i - 1][j];
        for (int k = 1; k <= count; k++)
        {
            int v = -1;
            if (k * weight <= j)
            {
                v = m[i - 1][j - k * weight] + k * value;
            }
            if (v > max)
            {
                max = v;
            }
        }

        m[i][j] = max;
    }
}

__efficient__ void knapsack_inner(int weight, int value, int count, int **m,
                                  int i, int w)
{
    knapsack_inner_loop(weight, value, count, m, i, w, 0, 8);
    knapsack_inner_loop(weight, value, count, m, i, w, 1, 8);
    knapsack_inner_loop(weight, value, count, m, i, w, 2, 8);
    knapsack_inner_loop(weight, value, count, m, i, w, 3, 8);
    knapsack_inner_loop(weight, value, count, m, i, w, 4, 8);
    knapsack_inner_loop(weight, value, count, m, i, w, 5, 8);
    knapsack_inner_loop(weight, value, count, m, i, w, 6, 8);
    knapsack_inner_loop(weight, value, count, m, i, w, 7, 8);
}

__efficient__ void knapsack_res(int *weights, int *values, int w, int *s,
                                int **m, int n)
{
    int i, j, k;
    j = w;
    for (i = n; i > 0; i--)
    {
        int v = m[i][j];
        int count = 0;
        for (k = 0; v != m[i - 1][j] + k * values[i - 1]; k++)
        {
            count++;
            j -= weights[i - 1];
        }

        s[i - 1] = count;
    }
}
#else
__efficient__ void knapsack_inner(int weight, int value, int count, int **m,
                                  int i, int w)
{
    int j, k, v;
    for (j = 0; j < w; j++)
    {
        m[i][j] = m[i - 1][j];
        for (k = 1; k <= count; k++)
        {
            if (k * weight > j)
            {
                break;
            }
            v = m[i - 1][j - k * weight] + k * value;
            if (v > m[i][j])
            {
                m[i][j] = v;
            }
        }
    }
}

__efficient__ void knapsack_res(int *weights, int *values, int w, int *s,
                                int **m, int n)
{
    int i, j, k;
    j = w;
    for (i = n; i > 0; i--)
    {
        int v = m[i][j];
        s[i - 1] = 0;
        for (k = 0; v != m[i - 1][j] + k * values[i - 1]; k++)
        {
            s[i - 1]++;
            j -= weights[i - 1];
        }
    }
}

#endif

void knapsack(int *weights, int *values, int *count, int w, int *mm, int **m,
              int *s, int n)
{
    int i;
    m[0] = mm;
    for (i = 1; i <= n; i++)
    {
        m[i] = &mm[i * (w + 1)];
        knapsack_inner(weights[i - 1], values[i - 1], count[i - 1], m, i,
                       w + 1);
    }
    knapsack_res(weights, values, w, s, m, n);
}
