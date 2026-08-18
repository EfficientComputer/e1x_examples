#ifdef __EFFCC__
#include <eff.h>
#endif

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MIN3(x, y, z) (MIN((x), MIN((y), (z))))

#ifndef EFF_BLD_HAND_OPTIMIZED

int dist[128][128];

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

#else

#define LEV_H 5

#define LEV_L2MAX 128
#define LEV_SLACK 32      // >= LEV_H - 1, at both ends
#define LEV_INF (1 << 20) // poison for a not-yet-valid diagonal
#define LEV_NOCHAR (-1)   // never equals an (unsigned char)

// Two explicit min2's, never MIN3 -- see the note above.
#define LEV_MIN2(a, b) ((a) < (b) ? (a) : (b))

// Load-bearing, not a hint: lane state must land in registers. Any dynamic
// index spills p/pp/s1c/c and costs ~4*LEV_H memory ops per step. Gate: no
// surviving `alloca` at -O3, constant getelementptr offsets.
#define LEV_UNROLL _Pragma("clang loop unroll(full)")

static int levRowBuf[LEV_SLACK + LEV_L2MAX + 1 + LEV_SLACK];

__attribute__((always_inline)) static void lev_band(const char *restrict s1,
                                                    const char *restrict s2,
                                                    int l2, int *restrict row,
                                                    int i0, int notFirst)
{
    int p[LEV_H];   // lane k's previous output, f[i0+k][j-k-1]
    int pp[LEV_H];  // lane k's output from two steps ago (lane k+1's diagonal)
    int s1c[LEV_H]; // band's s1 characters, loop-invariant
    int c[LEV_H];   // s2 shift register: c[k] == s2[j-k-1]
    int ptop;       // row[j-1] from the previous step

    LEV_UNROLL
    for (int k = 0; k < LEV_H; k++)
    {
        int rowIdx = i0 + k;
        // Virtual lanes (rowIdx <= 0) clone row 0, whose f is identically 0.
        p[k] = rowIdx > 0 ? rowIdx : 0; // f[i0+k][0]
        pp[k] = LEV_INF;
        s1c[k] = rowIdx > 0 ? (int)(unsigned char)s1[rowIdx - 1] : LEV_NOCHAR;
        c[k] = LEV_NOCHAR;
    }
    ptop = i0 > 1 ? i0 - 1 : 0; // f[i0-1][0]

    __effcc_ignore_memory_order
    {
        for (int j = 1; j < l2 + LEV_H; j++)
        {
            int v[LEV_H];

            int jj = j - 1 < l2 ? j - 1 : l2 - 1;
            LEV_UNROLL
            for (int k = LEV_H - 1; k > 0; k--)
                c[k] = c[k - 1];
            c[0] = (int)(unsigned char)s2[jj];

            int top0 = notFirst ? row[j] : 0;
            v[0] = LEV_MIN2(LEV_MIN2(top0 + 1, ptop - (s1c[0] == c[0])), p[0]);

            LEV_UNROLL
            for (int k = 1; k < LEV_H; k++)
                v[k] = LEV_MIN2(
                    LEV_MIN2(p[k - 1] + 1, pp[k - 1] - (s1c[k] == c[k])), p[k]);

            ptop = top0;

            LEV_UNROLL
            for (int k = 0; k < LEV_H; k++)
                pp[k] = (j > k) ? p[k] : LEV_INF;
            LEV_UNROLL
            for (int k = 0; k < LEV_H; k++)
                p[k] = v[k];

            row[j - LEV_H + 1] = v[LEV_H - 1];
        }
    }
}

__effcc_rip_exact void levenshtein_wave(const char *restrict s1,
                                        const char *restrict s2, int l2,
                                        int *restrict res, int nb, int r)
{
    int *restrict row = levRowBuf + LEV_SLACK;
    for (int b = 0; b < nb; b++)
        lev_band(s1, s2, l2, row, r - LEV_H + 1 + b * LEV_H, b);
    *res = row[l2] + l2; // undo the -j offset
}

void levenshtein(const char *restrict s1, int l1, const char *restrict s2,
                 int l2, int *restrict res)
{
    // Degenerate cases, kept off the fabric
    if (l1 <= 0)
    {
        *res = l2;
        return;
    }
    if (l2 <= 0)
    {
        *res = l1;
        return;
    }

    int nb = (l1 + LEV_H - 1) / LEV_H; // number of bands
    int r = l1 - (nb - 1) * LEV_H;     // rows in band 0, in [1, LEV_H]
    levenshtein_wave(s1, s2, l2, res, nb, r);
}

#endif
