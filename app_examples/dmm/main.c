#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "dmm.h"

#define MAT_REF_SIZE 128

int32_t _Alignas(16) mat_a[MAT_REF_SIZE][MAT_REF_SIZE];
int32_t _Alignas(16) mat_b[MAT_REF_SIZE][MAT_REF_SIZE];
int32_t _Alignas(16) mat_z[MAT_REF_SIZE][MAT_REF_SIZE];
int32_t _Alignas(16) mat_z_ref[MAT_REF_SIZE][MAT_REF_SIZE];

int randPrev = 42;
int pseudo_rand() {
    randPrev = (randPrev * 1000693) & 0x7FFF;
    return randPrev;
}

void generate_random_matrix_i32(int maxval, int minval, int n, int m, int32_t mat[restrict n][m])
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            mat[i][j] = (pseudo_rand() % (maxval - minval + 1)) + minval;
        }
    }
}

int main() {
    generate_random_matrix_i32(20, -20, MAT_REF_SIZE, MAT_REF_SIZE, mat_a);
    generate_random_matrix_i32(20, -20, MAT_REF_SIZE, MAT_REF_SIZE, mat_b);

    dmm_reference((int32_t *)mat_a, (int32_t *)mat_b, (int32_t *)mat_z_ref, MAT_REF_SIZE,
        MAT_REF_SIZE, MAT_REF_SIZE);

    dmm((int32_t *)mat_a, (int32_t *)mat_b, (int32_t *)mat_z, MAT_REF_SIZE,
        MAT_REF_SIZE, MAT_REF_SIZE);

    int32_t checksum = 0;
    for (uint32_t i = 0; i < MAT_REF_SIZE; i++) {
        for (uint32_t j = 0; j < MAT_REF_SIZE; j++) {
            if (mat_z[i][j] != mat_z_ref[i][j]) {
                printf("FAIL\r\n");
                return 0;
            }
        }
    }

    printf("PASS\r\n");
    return 0;
}
