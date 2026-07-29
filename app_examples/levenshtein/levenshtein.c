#ifdef __EFFCC__
#include <eff.h>
#endif

#define MIN(x, y) \
    ((x) < (y) ? (x) : (y)) // calculate minimum between two
                            // values

#define MIN3(x, y, z) \
    (MIN((x), MIN((y), (z)))) // calculate minimum between three
                              // values

int dist[128][128];

#ifdef EFF_BLD_HAND_OPTIMIZED
__attribute__((always_inline)) void levenshtein_inner(const char *s1, int l1,
                                                      const char *s2, int l2,
                                                      int *res, int startI)
{
    int pz = dist[startI - 1][0];
    int pp0 = 0x7FFFFFFE;
    int p0 = dist[startI][0];
    int pp1 = 0x7FFFFFFE;
    int p1 = dist[startI + 1][0];
    int pp2 = 0x7FFFFFFE;
    int p2 = dist[startI + 2][0];
    int pp3 = 0x7FFFFFFE;
    int p3 = dist[startI + 3][0];
    int p4 = dist[startI + 4][0];

    // 5 is the unroll amount
    for (int j = 1; j < l2 + 5; j++)
    {
        /*
         * Computes a diagonal of the matrix. Unrolled 5 times.
         * +-----
         * |    0 ->
         * |   1 ->
         * |  2 ->
         * | 3 ->
         * |4 ->
         *
         * For 1, the previous value of 0 is the top, its own previous value is
         * the left, and the one before the previous of 0 is the top left.
         */

        int a0 = j <= l2;
        int a1 = j > 1 && j <= l2 + 1 && startI + 1 <= l1;
        int a2 = j > 2 && j <= l2 + 2 && startI + 2 <= l1;
        int a3 = j > 3 && j <= l2 + 3 && startI + 3 <= l1;
        int a4 = j > 4 && startI + 4 <= l1;

        const char *s1p0 = &s1[startI - 1];
        const char *s1p1 = &s1[startI];
        const char *s1p2 = &s1[startI + 1];
        const char *s1p3 = &s1[startI + 2];
        const char *s1p4 = &s1[startI + 3];

        const char *s2p0 = &s2[j - 1];
        const char *s2p1 = &s2[j - 2];
        const char *s2p2 = &s2[j - 3];
        const char *s2p3 = &s2[j - 4];
        const char *s2p4 = &s2[j - 5];

        int t0 = a0 && *s1p0 == *s2p0 ? 0 : 1;
        int t1 = a1 && *s1p1 == *s2p1 ? 0 : 1;
        int t2 = a2 && *s1p2 == *s2p2 ? 0 : 1;
        int t3 = a3 && *s1p3 == *s2p3 ? 0 : 1;
        int t4 = a4 && *s1p4 == *s2p4 ? 0 : 1;

        int *vdp = &dist[startI - 1][j];
        int vz = j < l2 ? *vdp : 0x7FFFFFFE;
        int v0 = MIN3(vz + 1, pz + t0, p0 + 1);
        int v1 = MIN3(p1 + 1, p0 + 1, pp0 + t1);
        int v2 = MIN3(p2 + 1, p1 + 1, pp1 + t2);
        int v3 = MIN3(p3 + 1, p2 + 1, pp2 + t3);
        int v4 = MIN3(p4 + 1, p3 + 1, pp3 + t4);

        // updating the previous-previous values
        pp0 = p0;
        pp1 = p1;
        pp2 = p2;
        pp3 = p3;

        // updating the previous values
        pz = vz;

        p0 = v0;

        if (j > 1)
        {
            p1 = v1;
        }
        if (j > 2)
        {
            p2 = v2;
        }
        if (j > 3)
        {
            p3 = v3;
        }
        if (j > 4)
        {
            p4 = v4;
        }

        // storing results
        int *rdp0 = &dist[startI][j];
        int *rdp1 = &dist[startI + 1][j - 1];
        int *rdp2 = &dist[startI + 2][j - 2];
        int *rdp3 = &dist[startI + 3][j - 3];
        int *rdp4 = &dist[startI + 4][j - 4];

        if (a0)
            *rdp0 = v0;
        if (a1)
            *rdp1 = v1;
        if (a2)
            *rdp2 = v2;
        if (a3)
            *rdp3 = v3;
        if (a4)
            *rdp4 = v4;
    }
}

__efficient__ void levenshtein_inner_single(const char *s1, int l1,
                                            const char *s2, int l2, int *res,
                                            int startI)
{
    levenshtein_inner(s1, l1, s2, l2, res, startI);
}

__efficient__ void levenshtein_inner_full(const char *s1, int l1,
                                          const char *s2, int l2, int *res,
                                          int l1PlusOne)
{
    __effcc_region_exact for (int startI = 1; startI < l1PlusOne; startI += 5)
    {
        levenshtein_inner(s1, l1, s2, l2, res, startI);
    }
}

__efficient__ void levenshtein_init(const char *s1, int l1, const char *s2,
                                    int l2, int *res)
{
    for (int i = 0; i <= l1; i++)
    {
        dist[i][0] = i;
    }
    for (int j = 1; j <= l2; j++)
    {
        dist[0][j] = j;
    }
}

void levenshtein(const char *s1, int l1, const char *s2, int l2, int *res)
{
    if (l1 + l2 < 50)
    {
        // too short to overcome the overhead of running on the fabric
        for (int i = 0; i <= l1; i++)
        {
            dist[i][0] = i;
        }
        for (int j = 1; j <= l2; j++)
        {
            dist[0][j] = j;
        }
    }
    else
    {
        levenshtein_init(s1, l1, s2, l2, res);
    }

    for (int i = 1; i <= l1; i++)
    {
        levenshtein_inner_single(s1, l1, s2, l2, res, i);
    }

    *res = dist[l1][l2];
}
#else

__efficient__ void levenshtein(const char *s1, int l1, const char *s2, int l2,
                               int *res)
{
    for (int i = 0; i <= l1; i++)
    {
        dist[0][i] = i;
    }
    for (int j = 0; j <= l2; j++)
    {
        dist[j][0] = j;
    }
    for (int j = 1; j <= l1; j++)
    {
        for (int i = 1; i <= l2; i++)
        {
            int track;
            if (s2[i - 1] == s1[j - 1])
            {
                track = 0;
            }
            else
            {
                track = 1;
            }
            int t = MIN((dist[i - 1][j] + 1), (dist[i][j - 1] + 1));
            dist[i][j] = MIN(t, (dist[i - 1][j - 1] + track));
        }
    }
    *res = dist[l2][l1];
}

#endif
