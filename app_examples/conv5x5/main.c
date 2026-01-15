#include <eff/profile.h>
#include <stdlib.h>
#include "conv5x5.h"

#define N 64

#define N_PAD (N + 4)
#define RANDOMIZE_FILTER 1
#define RANGE 10

data_t filter_global[5 * 5] = {
    DATA_INIT(0),  DATA_INIT(-1), DATA_INIT(2), DATA_INIT(-2), DATA_INIT(0),
    DATA_INIT(-1), DATA_INIT(-2), DATA_INIT(3), DATA_INIT(-2), DATA_INIT(-1),
    DATA_INIT(0),  DATA_INIT(-3), DATA_INIT(5), DATA_INIT(-3), DATA_INIT(0),
    DATA_INIT(-2), DATA_INIT(-3), DATA_INIT(5), DATA_INIT(2),  DATA_INIT(1),
    DATA_INIT(-1), DATA_INIT(0),  DATA_INIT(3), DATA_INIT(1),  DATA_INIT(-1)};

data_t in[N_PAD * N_PAD];
data_t ref_out[N_PAD * N_PAD];
data_t test_out[N * N];

void dconv_ref(const data_t *in, int n, const data_t *filter, data_t *out) {
    const int d = 5;
    for (int i = 0; i <= n - d; i++) {
        for (int j = 0; j <= n - d; j++) {
            data_t sum = DATA_INIT(0);

            for (int fi = 0; fi < d; fi++) {
                for (int fj = 0; fj < d; fj++) {
                    int col = i + fi;
                    int row = j + fj;

                    sum += filter[fi + fj * d] * in[col + row * n];
                }
            }

            out[i + j * n] = sum;
        }
    }
}

void print_matrix(const data_t *mat, int n, int stride) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%08x ", mat[i * stride + j]);
        }
        printf("\n");
    }
}

int randPrev = 42;
int pseudo_rand() {
    randPrev = (randPrev * 1000693) & 0x7FFF;
    return randPrev;
}

int main() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            ref_out[i + j * N] = DATA_INIT(0);
            test_out[i + j * N] = DATA_INIT(0);
            in[i + j * N] = DATA_INIT(pseudo_rand()) % RANGE;
        }
    }

    for (int i = 0; i < 5 * 5; i++) {
        int f = pseudo_rand() % RANGE - (RANGE - 1) / 2;
        filter_global[i] = DATA_INIT(f);
    }

    START_PROFILE_REGION("kernel");
    for (int iter = 0; iter < 100; iter++)
        dconv5x5(in, N, N, filter_global, test_out);
    END_PROFILE_REGION();

    dconv_ref(in, N, filter_global, ref_out);

    for (int i = 0; i < N - 4; i++) {
        for (int j = 0; j < N - 4; j++) {
            if ((int)ref_out[i + j * N] != (int)test_out[i + j * N]) {
                printf("[conv5x5] FAIL i: %d j: %d - %d != %d\n", i, j,
                       (int)ref_out[i + j * N], (int)test_out[i + j * N]);
                return 1;
            }
        }
    }

    printf("[conv5x5] PASS\n");
    return 0;
}
