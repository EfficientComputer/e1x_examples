#include <stdio.h>
#include "qr_decomp.h"
#include "qr.h"
#include <stdint.h>
#include <stdlib.h>

matrix *a;
matrix *q;
matrix *r;
matrix *v;

int pseudo_rand(int prev) { return (prev * 1000693) & 0x7FFF; }

void setupQRDecomp() {
    a = makeMatrix(INPUT_SIZE, INPUT_SIZE);
    q = makeMatrix(INPUT_SIZE, INPUT_SIZE);
    r = makeMatrix(INPUT_SIZE, INPUT_SIZE);

    int prev = pseudo_rand(12 /* = seed */);
    for (int i = 0; i < INPUT_SIZE; ++i) {
        for (int j = 0; j < INPUT_SIZE; ++j) {
            prev = pseudo_rand(prev);
            a->data[i * INPUT_SIZE + j] = FIXED_POINT_CAST_TO(prev % 15 - 7);
        }
    }

    v = NULL;
}

void QRDecomp() {
    if (v != NULL) {
        freeMatrix(v);
    }
    v = copyMatrix(a);
    gram_schmidt(v, &q, &r);
}

void checkQRDecomp() {
    int expected_checksum = EXPECTED_CHECKSUM;

    int checksum = 0;
    for (int i = 0; i < INPUT_SIZE; ++i) {
        for (int j = 0; j < INPUT_SIZE; ++j) {
            checksum ^=
                q->data[i * INPUT_SIZE + j] + r->data[i * INPUT_SIZE + j];
        }
    }

    if (checksum == expected_checksum) {
        printf("[qr_decomp] PASS\n");
    } else {
        printf("[qr_decomp] FAIL. Got checksum %d, expected %d\n", checksum,
               expected_checksum);
    }

    freeMatrix(a);
    freeMatrix(v);
    freeMatrix(q);
    freeMatrix(r);
}