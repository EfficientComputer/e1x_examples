#include <stdint.h>
#ifdef __EFFCC__
#include <effintrinsic.h>
#endif

// 4-row unroll with vdot_s16: inline int32->int16 cast, step j by 2.
// Assumes n divisible by 4, m divisible by 2.
__efficient__ void dmv(int32_t *a, int32_t *b, int32_t *z, uint32_t n,
                       uint32_t m) {
    for (uint32_t i = 0; i < n; i += 4) {
        int32_t w0 = 0, w1 = 0, w2 = 0, w3 = 0;
        const int32_t *a0 = a + (i + 0) * m;
        const int32_t *a1 = a + (i + 1) * m;
        const int32_t *a2 = a + (i + 2) * m;
        const int32_t *a3 = a + (i + 3) * m;
        __effcc_ignore_memory_order for (uint32_t j = 0; j < m; j += 2) {
            _v2s16_t bv = {(int16_t)b[j], (int16_t)b[j + 1]};
            _v2s16_t a0v = {(int16_t)a0[j], (int16_t)a0[j + 1]};
            _v2s16_t a1v = {(int16_t)a1[j], (int16_t)a1[j + 1]};
            _v2s16_t a2v = {(int16_t)a2[j], (int16_t)a2[j + 1]};
            _v2s16_t a3v = {(int16_t)a3[j], (int16_t)a3[j + 1]};
            w0 += __builtin_effcc_vdot_s16(a0v, bv);
            w1 += __builtin_effcc_vdot_s16(a1v, bv);
            w2 += __builtin_effcc_vdot_s16(a2v, bv);
            w3 += __builtin_effcc_vdot_s16(a3v, bv);
        }
        z[i + 0] = w0;
        z[i + 1] = w1;
        z[i + 2] = w2;
        z[i + 3] = w3;
    }
}
