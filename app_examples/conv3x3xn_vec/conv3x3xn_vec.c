#ifdef __EFFCC__
#include <effintrinsic.h>
#endif

#include "conv3x3xn_vec.h"

// vdot_s16 vectorized multi-channel 3x3 conv.
// k (channel) is the outer loop so per-channel weight packing is done once
// per channel, not once per output pixel. Output is zero-initialised then
// accumulated with +=. No memory-ordering relaxation: sequential k loop
// serialises cross-channel writes to the same output location.
__efficient__ void dconv_3x3xn_internal(const data_t *in, const data_t *kernels,
                                        int nchannels, data_t *out)
{
    for (int idx = 0; idx < M * OUTSTRIDE; idx++)
        out[idx] = 0;

    for (int k = 0; k < nchannels; ++k)
    {
        const data_t *ks = kernels + KERNEL_DIM * KERNEL_DIM * k;

        // Pack channel k's 9 taps once — hoisted above the i/j pixel loops.
        _v2s16_t kv0 = {(int16_t)ks[0], (int16_t)ks[1]};
        _v2s16_t kv1 = {(int16_t)ks[2], (int16_t)ks[3]};
        _v2s16_t kv2 = {(int16_t)ks[4], (int16_t)ks[5]};
        _v2s16_t kv3 = {(int16_t)ks[6], (int16_t)ks[7]};
        int16_t k8 = (int16_t)ks[8];

        for (int i = 0; i < M - KERNEL_DIM + 1; i++)
        {
            for (int j = 0; j < N - KERNEL_DIM + 1; j++)
            {
                _v2s16_t iv0 = {(int16_t)in[(i + 0) * INSTRIDE + (j + 0)],
                                (int16_t)in[(i + 0) * INSTRIDE + (j + 1)]};
                _v2s16_t iv1 = {(int16_t)in[(i + 0) * INSTRIDE + (j + 2)],
                                (int16_t)in[(i + 1) * INSTRIDE + (j + 0)]};
                _v2s16_t iv2 = {(int16_t)in[(i + 1) * INSTRIDE + (j + 1)],
                                (int16_t)in[(i + 1) * INSTRIDE + (j + 2)]};
                _v2s16_t iv3 = {(int16_t)in[(i + 2) * INSTRIDE + (j + 0)],
                                (int16_t)in[(i + 2) * INSTRIDE + (j + 1)]};

                out[i * OUTSTRIDE + j] +=
                    __builtin_effcc_vdot_s16(kv0, iv0) +
                    __builtin_effcc_vdot_s16(kv1, iv1) +
                    __builtin_effcc_vdot_s16(kv2, iv2) +
                    __builtin_effcc_vdot_s16(kv3, iv3) +
                    (data_t)k8 *
                        (data_t)(int16_t)in[(i + 2) * INSTRIDE + (j + 2)];
            }
        }
    }
}

void dconv_3x3xn(const data_t *in, int m, int n, int instride, int outstride,
                 const data_t *kernels, int nchannels, data_t *out)
{
    (void)m;
    (void)n;
    (void)instride;
    (void)outstride;
    dconv_3x3xn_internal(in, kernels, nchannels, out);
}
