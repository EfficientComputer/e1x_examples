#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#ifndef EFFX_INPUT_SIZE
#define EFFX_INPUT_SIZE 2
#endif

#if EFFX_INPUT_SIZE == 2
#define MAT_REF_SIZE 256
#elif EFFX_INPUT_SIZE == 1
#define MAT_REF_SIZE 128
#elif EFFX_INPUT_SIZE == 0
#define MAT_REF_SIZE 64
#endif

#if MAT_REF_SIZE % 4 != 0
#error "MAT_REF_SIZE must be divisible by 4"
#endif

int32_t mat_a[MAT_REF_SIZE][MAT_REF_SIZE] __attribute__((aligned(16)));
int32_t mat_b[MAT_REF_SIZE][MAT_REF_SIZE] __attribute__((aligned(16)));
int32_t mat_z[MAT_REF_SIZE][MAT_REF_SIZE] __attribute__((aligned(16)));
int32_t mat_z_ref[MAT_REF_SIZE][MAT_REF_SIZE] __attribute__((aligned(16)));

void dmm(int32_t *a,  // n x m
         int32_t *b,  // m x o
         int32_t *z,  // n x o
         uint32_t n, uint32_t m, uint32_t o);

void dmm_reference(
         int32_t *a,  // n x m
         int32_t *b,  // m x o
         int32_t *z,  // n x o
         uint32_t n, uint32_t m, uint32_t o);

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
                printf("FAIL\n");
                return 0;
            }
        }
    }

    printf("PASS\n");
    return 0;
}
