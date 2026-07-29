#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"

void assert(int assertion, char *message)
{
    if (assertion == 0)
    {
        printf("%s\n", message);
        return; // exit(1);
    }
}

matrix *makeMatrix(int width, int height)
{
    matrix *out;
    assert(width > 0 && height > 0, "New matrix must be at least a 1 by 1");
    out = (matrix *)malloc(sizeof(matrix));

    assert(out != NULL, "Out of memory.");

    out->width = width;
    out->height = height;
    out->data = (int *)malloc(sizeof(int) * width * height);

    assert(out->data != NULL, "Out of memory.");

    memset(out->data, 0, width * height * sizeof(int));

    return out;
}

void freeMatrix(matrix *m)
{
    if (m != NULL)
    {
        if (m->data != NULL)
        {
            free(m->data);
            m->data = NULL;
        }

        free(m);
        m = NULL;
    }
    return;
}

matrix *scaleMatrix(matrix *m, int value)
{
    int i, elements = m->width * m->height;
    matrix *out = makeMatrix(m->width, m->height);
    int *ptrM = m->data;
    int *ptrOut = out->data;

    for (i = 0; i < elements; i++)
    {
        *(ptrOut++) = *(ptrM++) * value;
    }

    return out;
}

__attribute__((always_inline)) matrix *copyMatrix(matrix *m)
{
    return scaleMatrix(m, 1);
}

void printMatrix(matrix *m)
{
    int i, j;
    int *ptr = m->data;
    printf("%d %d\n", m->width, m->height);
    for (i = 0; i < m->height; i++)
    {
        for (j = 0; j < m->width; j++)
        {
            printf(" %d", *(ptr++));
        }
        printf("\n");
    }
    return;
}
