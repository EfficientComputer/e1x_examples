#ifdef __EFFCC__
#include <eff.h>
#endif

#include "qr.h"

#define MAT_AT(matrix, row, col) \
    ((matrix)->data + (matrix)->width * (row) + (col))

__attribute__((always_inline)) int32_t isqrt(int32_t x)
{
    int32_t q = 1, r = 0;
    while (q <= x)
    {
        q <<= 2;
    }
    while (q > 1)
    {
        int32_t t;
        q >>= 2;
        t = x - r - q;
        r >>= 1;
        if (t >= 0)
        {
            x = t;
            r += q;
        }
    }
    return r;
}

#ifdef EFF_BLD_HAND_OPTIMIZED
__efficient__ void update_q_and_r(int *restrict vdata, int *restrict qdata,
                                  int *restrict rdata, int width, int height,
                                  int kp1)
{
    __effcc_ignore_memory_order
    {
        for (int j = kp1; j < width; j++)
        {
            int sum = 0;
            for (int i = 0; i < height; i++)
            {
                sum += FIXED_POINT_MUL_RESCALE(vdata[i * width + j] *
                                               qdata[i * width + kp1 - 1]);
            }

            for (int i = 0; i < height; i++)
            {
                vdata[i * width + j] -=
                    FIXED_POINT_MUL_RESCALE(qdata[i * width + kp1 - 1] * sum);
            }

            rdata[(kp1 - 1) * width + j] = sum;
        }
    }
}

void gram_schmidt(matrix *v, matrix **q, matrix **r)
{
    int i, j, k;
    int l2norm;
    int sum;
    int *vPtr;
    int *qPtr;
    int *rPtr;

    // For each column in A (now called V)
    for (k = 0; k < v->width; k++)
    {
        vPtr = MAT_AT(v, 0, k);
        // Step 1: Get the L2-Norm of column k
        l2norm = 0;
        for (i = 0; i < v->height; i++)
        {
            l2norm += FIXED_POINT_MUL_RESCALE(*vPtr * *vPtr);
            vPtr += v->width;
        }

        // Take the square root interpreting the fixed-point value as
        // a regular 32-bit int. Then divide by the sqare root of the FP
        // scaling factor.
        l2norm = FIXED_POINT_DIVIDEND(isqrt(l2norm)) >> (FIXED_POINT_BITS / 2);

        // Store this value in R(k,k)
        // The nice thing about the rPtr variable is that
        // it only has to be readjusted at the beginning of the
        // first 'for' loop. After each use, we just increment by 1.
        rPtr = MAT_AT(*r, k, k);
        *rPtr = l2norm;
        rPtr++;

        // Step 2: Normalize A's column k and store it in Q.
        vPtr = MAT_AT(v, 0, k);
        qPtr = MAT_AT(*q, 0, k);
        for (i = 0; i < v->height; i++)
        {
            *qPtr = (l2norm == 0 ? 0 : FIXED_POINT_DIVIDEND(*vPtr) / l2norm);
            vPtr += v->width;
            qPtr += (*q)->width;
        }

        int *vdata = v->data;
        int *qdata = (*q)->data;
        int *rdata = (*r)->data;
        int width = v->width;
        int height = v->height;
        update_q_and_r(vdata, qdata, rdata, width, height, k + 1);
    }
}
#else
__efficient__ void gram_schmidt(matrix *v, matrix **q, matrix **r)
{
    int i, j, k;
    int l2norm;
    int sum;
    int *vPtr;
    int *qPtr;
    int *rPtr;

    // For each column in A (now called V)
    for (k = 0; k < v->width; k++)
    {
        vPtr = MAT_AT(v, 0, k);
        // Step 1: Get the L2-Norm of column k
        l2norm = 0;
        for (i = 0; i < v->height; i++)
        {
            l2norm += FIXED_POINT_MUL_RESCALE(*vPtr * *vPtr);
            vPtr += v->width;
        }

        // Take the square root interpreting the fixed-point value as
        // a regular 32-bit int. Then divide by the sqare root of the FP
        // scaling factor.
        l2norm = FIXED_POINT_DIVIDEND(isqrt(l2norm)) >> (FIXED_POINT_BITS / 2);

        // Store this value in R(k,k)
        // The nice thing about the rPtr variable is that
        // it only has to be readjusted at the beginning of the
        // first 'for' loop. After each use, we just increment by 1.
        rPtr = MAT_AT(*r, k, k);
        *rPtr = l2norm;
        rPtr++;

        // Step 2: Normalize A's column k and store it in Q.
        vPtr = MAT_AT(v, 0, k);
        qPtr = MAT_AT(*q, 0, k);
        for (i = 0; i < v->height; i++)
        {
            *qPtr = (l2norm == 0 ? 0 : FIXED_POINT_DIVIDEND(*vPtr) / l2norm);
            vPtr += v->width;
            qPtr += (*q)->width;
        }

        // Step 3: 2 parts. For each column after K, do the following:
        for (j = k + 1; j < v->width; j++)
        {
            // Step 3a: Dot Product Q's column K with A's column J,
            // storing the result of the Dot Product at R(k,j)
            qPtr = MAT_AT(*q, 0, k);
            vPtr = MAT_AT(v, 0, j);
            sum = 0;
            for (i = 0; i < v->height; i++)
            {
                sum += FIXED_POINT_MUL_RESCALE(*vPtr * *qPtr);
                vPtr += v->width;
                qPtr += (*q)->width;
            }
            *rPtr = sum;
            rPtr++;

            // Step 3b: Multiply Q's column K with sum
            // (which is stored at *rPtr). Then take A's column J
            // and subtract from it Q's K * sum. Take this
            // result and store it back into A's column J.
            vPtr = MAT_AT(v, 0, j);
            qPtr = MAT_AT(*q, 0, k);
            for (i = 0; i < v->height; i++)
            {
                *vPtr = *vPtr - FIXED_POINT_MUL_RESCALE(*qPtr * sum);
                vPtr += v->width;
                qPtr += (*q)->width;
            }
        }
    }
}
#endif
