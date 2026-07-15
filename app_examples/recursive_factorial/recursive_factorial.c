__efficient__ void factorial(int n, int *y)
{
    if (n != 0)
    {
        *y = *y * n;
        factorial(n - 1, y);
    }
}
