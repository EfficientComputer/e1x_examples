// <<AUTOBENCH>> efficient_e1 efficient_e1x
// clang-format off
#include "stopprop.h" // <<AUTOBENCH-INCLUDE>> ../../../include/stopprop.h
// clang-format on
#define DEFINE_TWIDDLES
#include "fft4k.h"

#define FIXED_ROUND(x) \
    (eff_fft_scalar)(((x) + (1 << (FRACBITS - 1))) >> FRACBITS)

#define SAMP_MAX 32767
#define FRACBITS 15

// Naive implementation macros
#define SAMPPROD int32_t
#define smul(a, b) ((SAMPPROD)(a) * (b))
#define sround(x) (eff_fft_scalar)(((x) + (1 << (FRACBITS - 1))) >> FRACBITS)

#define S_MUL(a, b) sround(smul(a, b))

#define C_MUL(m, a, b)                                           \
    do {                                                         \
        (m).r = sround(smul((a).r, (b).r) - smul((a).i, (b).i)); \
        (m).i = sround(smul((a).r, (b).i) + smul((a).i, (b).r)); \
    } while (0)

#define DIVSCALAR(x, k) (x) = sround(smul(x, SAMP_MAX / k))

#define C_FIXDIV(c, div)       \
    do {                       \
        DIVSCALAR((c).r, div); \
        DIVSCALAR((c).i, div); \
    } while (0)

#define C_MULBYSCALAR(c, s)             \
    do {                                \
        (c).r = sround(smul((c).r, s)); \
        (c).i = sround(smul((c).i, s)); \
    } while (0)

#ifndef CHECK_OVERFLOW_OP
#define CHECK_OVERFLOW_OP(a, op, b) /* noop */
#endif

#define C_ADD(res, a, b)                   \
    do {                                   \
        CHECK_OVERFLOW_OP((a).r, +, (b).r) \
        CHECK_OVERFLOW_OP((a).i, +, (b).i) \
        (res).r = (a).r + (b).r;           \
        (res).i = (a).i + (b).i;           \
    } while (0)
#define C_SUB(res, a, b)                   \
    do {                                   \
        CHECK_OVERFLOW_OP((a).r, -, (b).r) \
        CHECK_OVERFLOW_OP((a).i, -, (b).i) \
        (res).r = (a).r - (b).r;           \
        (res).i = (a).i - (b).i;           \
    } while (0)
#define C_ADDTO(res, a)                      \
    do {                                     \
        CHECK_OVERFLOW_OP((res).r, +, (a).r) \
        CHECK_OVERFLOW_OP((res).i, +, (a).i) \
        (res).r += (a).r;                    \
        (res).i += (a).i;                    \
    } while (0)

#define C_SUBFROM(res, a)                    \
    do {                                     \
        CHECK_OVERFLOW_OP((res).r, -, (a).r) \
        CHECK_OVERFLOW_OP((res).i, -, (a).i) \
        (res).r -= (a).r;                    \
        (res).i -= (a).i;                    \
    } while (0)

#ifndef EFF_ARCH_E0

// FFT init is a NOP for Efficient devices; only necessary for ARM FFT
void fft_init() { return; }

#ifdef EFF_BLD_HAND_OPTIMIZED
__efficient__ void eff_fft_init_dst(eff_fft_cpx* dst, eff_fft_cpx* src,
                                    int size) {
    __effcc_ignore_memory_order {
        // Computing log2(size)
        uint32_t log2Size = -1;
        uint32_t sizeShift = size;
        while (sizeShift) {
            sizeShift >>= 1;
            log2Size++;
        }

        int highShiftAmt = 32 - log2Size + 1;
        int lowShiftAmt = 32 - log2Size - 1;

        __effcc_parallel(16) {
            for (int i = 0; i < size; i++) {
                // Reversing bits of dstIdx (i)
                uint32_t reversed = __builtin_elementwise_bitreverse(i);

                // Swapping each pair of bits
                uint32_t reversedA = reversed & 0xAAAAAAAA;
                uint32_t reversedB = reversed & 0x55555555;

                reversedA = reversedA >> highShiftAmt;
                reversedB = reversedB >> lowShiftAmt;

                uint32_t srcIdx = reversedA | reversedB;

                dst[i] = src[srcIdx];
            }
        }
    }
}

__efficient__ void eff_fft_run_layer(eff_fft_cpx* data, int twiddleStart,
                                     int idxStride, int scheduleLen) {
    __effcc_ignore_memory_order {
        for (int i = 0; i < scheduleLen; i++) {
            uint32_t* Fout = (uint32_t*)(data + (i * idxStride));
            uint32_t m = idxStride / 4;

            uint32_t *tw1, *tw2, *tw3;
            uint32_t m2 = m + m;
            uint32_t m3 = m2 + m;

            tw3 = tw2 = tw1 = (uint32_t*)(twiddles + twiddleStart * 2);

            for (int32_t k = m; k > 0; k--) {
                // C_FIXDIV(*Fout,4);
                uint32_t Fout0 = Fout[0];
                int16_t Fout0r = (int16_t)(Fout0 & 0xFFFF);
                int16_t Fout0i = (int16_t)(Fout0 >> 16);
                int32_t fOut0r =
                    FIXED_ROUND((int32_t)(Fout0r * (SAMP_MAX / 4)));
                int32_t fOut0i =
                    FIXED_ROUND((int32_t)(Fout0i * (SAMP_MAX / 4)));
                // C_FIXDIV(Fout[m],4);
                uint32_t Fout1 = Fout[m];
                int16_t Fout1r = (int16_t)(Fout1 & 0xFFFF);
                int16_t Fout1i = (int16_t)(Fout1 >> 16);
                int32_t fOut1r =
                    FIXED_ROUND((int32_t)(Fout1r * (SAMP_MAX / 4)));
                int32_t fOut1i =
                    FIXED_ROUND((int32_t)(Fout1i * (SAMP_MAX / 4)));
                // C_FIXDIV(Fout[m2],4);
                uint32_t Fout2 = Fout[m2];
                int16_t Fout2r = (int16_t)(Fout2 & 0xFFFF);
                int16_t Fout2i = (int16_t)(Fout2 >> 16);
                int32_t fOut2r =
                    FIXED_ROUND((int32_t)(Fout2r * (SAMP_MAX / 4)));
                int32_t fOut2i =
                    FIXED_ROUND((int32_t)(Fout2i * (SAMP_MAX / 4)));
                // C_FIXDIV(Fout[m3],4);
                uint32_t Fout3 = Fout[m3];
                int16_t Fout3r = (int16_t)(Fout3 & 0xFFFF);
                int16_t Fout3i = (int16_t)(Fout3 >> 16);
                int32_t fOut3r =
                    FIXED_ROUND((int32_t)(Fout3r * (SAMP_MAX / 4)));
                int32_t fOut3i =
                    FIXED_ROUND((int32_t)(Fout3i * (SAMP_MAX / 4)));

                // C_MUL(scratch[0],Fout[m] , *tw1 );
                uint32_t twiddle1 = *tw1;
                int16_t twiddle1r = (int16_t)(twiddle1 & 0xFFFF);
                int16_t twiddle1i = (int16_t)(twiddle1 >> 16);
                int32_t scratch0r = FIXED_ROUND(
                    (int32_t)(fOut1r * twiddle1r - fOut1i * twiddle1i));
                int32_t scratch0i = FIXED_ROUND(
                    (int32_t)(fOut1r * twiddle1i + fOut1i * twiddle1r));

                // C_MUL(scratch[1],Fout[m2] , *tw2 );
                uint32_t twiddle2 = *tw2;
                int16_t twiddle2r = (int16_t)(twiddle2 & 0xFFFF);
                int16_t twiddle2i = (int16_t)(twiddle2 >> 16);
                int32_t scratch1r = FIXED_ROUND(
                    (int32_t)(fOut2r * twiddle2r - fOut2i * twiddle2i));
                int32_t scratch1i = FIXED_ROUND(
                    (int32_t)(fOut2r * twiddle2i + fOut2i * twiddle2r));

                // C_MUL(scratch[2],Fout[m3] , *tw3 );
                uint32_t twiddle3 = *tw3;
                int16_t twiddle3r = (int16_t)(twiddle3 & 0xFFFF);
                int16_t twiddle3i = (int16_t)(twiddle3 >> 16);
                int32_t scratch2r = FIXED_ROUND(
                    (int32_t)(fOut3r * twiddle3r - fOut3i * twiddle3i));
                int32_t scratch2i = FIXED_ROUND(
                    (int32_t)(fOut3r * twiddle3i + fOut3i * twiddle3r));

                // C_SUB( scratch[5] , *Fout, scratch[1] );
                int32_t scratch5r = fOut0r - scratch1r;
                int32_t scratch5i = fOut0i - scratch1i;

                // C_ADDTO(*Fout, scratch[1]);
                fOut0r = fOut0r + scratch1r;
                fOut0i = fOut0i + scratch1i;

                // C_ADD( scratch[3] , scratch[0] , scratch[2] );
                int32_t scratch3r = scratch0r + scratch2r;
                int32_t scratch3i = scratch0i + scratch2i;

                // C_SUB( scratch[4] , scratch[0] , scratch[2] );
                int32_t scratch4r = scratch0r - scratch2r;
                int32_t scratch4i = scratch0i - scratch2i;

                // C_SUB( Fout[m2], *Fout, scratch[3] );
                Fout[m2] = (uint32_t)((uint16_t)(fOut0r - scratch3r)) |
                           (((uint16_t)(fOut0i - scratch3i)) << 16);

                tw1 += 1;
                tw2 += 2;
                tw3 += 3;
                // C_ADDTO( *Fout , scratch[3] );
                Fout[0] = (uint32_t)((uint16_t)(fOut0r + scratch3r)) |
                          (((uint16_t)(fOut0i + scratch3i)) << 16);

#if FFT_IS_INVERSE == 1
                Fout[m] = (uint32_t)((uint16_t)(scratch5r - scratch4i)) |
                          (((uint16_t)(scratch5i + scratch4r)) << 16);
                Fout[m3] = (uint32_t)((uint16_t)(scratch5r + scratch4i)) |
                           (((uint16_t)(scratch5i - scratch4r)) << 16);
#else
                Fout[m] = (uint32_t)((uint16_t)(scratch5r + scratch4i)) |
                          (((uint16_t)(scratch5i - scratch4r)) << 16);
                Fout[m3] = (uint32_t)((uint16_t)(scratch5r - scratch4i)) |
                           (((uint16_t)(scratch5i + scratch4r)) << 16);
#endif
                ++Fout;
            }
        }
    }
}

#else  // EFF_BLD_HAND_OPTIMIZED

__efficient__ void eff_fft_init_dst(eff_fft_cpx* dst, eff_fft_cpx* src,
                                    int size) {
    // Computing log2(size)
    uint32_t log2Size = -1;
    uint32_t sizeShift = size;
    while (sizeShift) {
        sizeShift >>= 1;
        log2Size++;
    }

    int highShiftAmt = 32 - log2Size + 1;
    int lowShiftAmt = 32 - log2Size - 1;

    for (int i = 0; i < size; i++) {
        // Reversing bits of dstIdx (i)
        uint32_t reversed = __builtin_elementwise_bitreverse(i);

        // Swapping each pair of bits
        uint32_t reversedA = reversed & 0xAAAAAAAA;
        uint32_t reversedB = reversed & 0x55555555;

        reversedA = reversedA >> highShiftAmt;
        reversedB = reversedB >> lowShiftAmt;

        uint32_t srcIdx = reversedA | reversedB;

        dst[i] = src[srcIdx];
    }
}

__efficient__ void eff_fft_run_layer(eff_fft_cpx* data, int twiddleStart,
                                     int idxStride, int scheduleLen) {
    for (int i = 0; i < scheduleLen; i++) {
        eff_fft_cpx* Fout = (eff_fft_cpx*)(data + (i * idxStride));
        uint32_t m = idxStride / 4;

        eff_fft_cpx *tw1, *tw2, *tw3;
        uint32_t m2 = m + m;
        uint32_t m3 = m2 + m;

        tw3 = tw2 = tw1 = (eff_fft_cpx*)(twiddles + twiddleStart * 2);

        eff_fft_cpx scratch[6];

        for (int32_t k = m; k > 0; k--) {
            C_FIXDIV(*Fout, 4);
            C_FIXDIV(Fout[m], 4);
            C_FIXDIV(Fout[m2], 4);
            C_FIXDIV(Fout[m3], 4);

            C_MUL(scratch[0], Fout[m], *tw1);
            C_MUL(scratch[1], Fout[m2], *tw2);
            C_MUL(scratch[2], Fout[m3], *tw3);

            C_SUB(scratch[5], *Fout, scratch[1]);

            C_ADDTO(*Fout, scratch[1]);

            C_ADD(scratch[3], scratch[0], scratch[2]);

            C_SUB(scratch[4], scratch[0], scratch[2]);

            C_SUB(Fout[m2], *Fout, scratch[3]);

            tw1 += 1;
            tw2 += 2;
            tw3 += 3;

            C_ADDTO(*Fout, scratch[3]);

#if FFT_IS_INVERSE == 1
            Fout[m].r = scratch[5].r - scratch[4].i;
            Fout[m].i = scratch[5].i + scratch[4].r;
            Fout[m3].r = scratch[5].r + scratch[4].i;
            Fout[m3].i = scratch[5].i - scratch[4].r;
#else
            Fout[m].r = scratch[5].r + scratch[4].i;
            Fout[m].i = scratch[5].i - scratch[4].r;
            Fout[m3].r = scratch[5].r - scratch[4].i;
            Fout[m3].i = scratch[5].i + scratch[4].r;
#endif
            ++Fout;
        }
    }
}

#endif  // EFF_BLD_HAND_OPTIMIZED

void eff_fft4(eff_fft_cpx* src, eff_fft_cpx* dst) {
    eff_fft_init_dst(dst, src, FFT_SIZE);

    int i = 0;
    int maxStride = FFT_SIZE;
    for (int m = 1; m < maxStride; m *= 4) {
        int scheduleLen = maxStride / (m * 4);
        eff_fft_run_layer(dst, twiddleSchedule[i], m * 4, scheduleLen);
        i++;
    }
}

#else

__efficient__ void eff_fft_init_dst(eff_fft_cpx* dst, eff_fft_cpx* src,
                                    int size, int shiftAmt) {
    __effcc_parallel(6) {
        for (int i = 0; i < size; i++) {
            // Reversing bits of dstIdx (i)
            uint32_t reversed = __builtin_elementwise_bitreverse(i);

            uint32_t srcIdx = reversed >> shiftAmt;
            dst[i] = src[srcIdx];
        }
    }
}

__efficient__ void eff_fft_run_layer(eff_fft_cpx* data, int twiddleStart,
                                     int idxStride, int scheduleLen) {
    for (int i = 0; i < scheduleLen; i++) {
        uint32_t* Fout = (uint32_t*)(data + (i * idxStride));
        uint32_t m = idxStride >> 1;

        uint32_t* tw1;

        tw1 = (uint32_t*)(twiddles + twiddleStart * 2);

        for (int32_t k = m; k > 0; k--) {
            // C_FIXDIV(*Fout, 2);
            uint32_t Fout0 = Fout[0];
            int16_t Fout0r = (int16_t)(Fout0 & 0xFFFF);
            int16_t Fout0i = (int16_t)(Fout0 >> 16);
            int32_t fOut0r = FIXED_ROUND((int32_t)(Fout0r * (SAMP_MAX / 2)));
            int32_t fOut0i = FIXED_ROUND((int32_t)(Fout0i * (SAMP_MAX / 2)));
            // C_FIXDIV(*Fout[m], 2);
            uint32_t Fout1 = Fout[m];
            int16_t Fout1r = (int16_t)(Fout1 & 0xFFFF);
            int16_t Fout1i = (int16_t)(Fout1 >> 16);
            int32_t fOut1r = FIXED_ROUND((int32_t)(Fout1r * (SAMP_MAX / 2)));
            int32_t fOut1i = FIXED_ROUND((int32_t)(Fout1i * (SAMP_MAX / 2)));

            // C_MUL(t, *Fout2[m], *tw1);
            uint32_t twiddle1 = *tw1;
            int16_t twiddle1r = (int16_t)(twiddle1 & 0xFFFF);
            int16_t twiddle1i = (int16_t)(twiddle1 >> 16);
            int32_t scratch0r =
                FIXED_ROUND((int32_t)(fOut1r * twiddle1r - fOut1i * twiddle1i));
            int32_t scratch0i =
                FIXED_ROUND((int32_t)(fOut1r * twiddle1i + fOut1i * twiddle1r));

            tw1 += 1;
            // C_SUB(*Fout2[m], *Fout, t);
            Fout[m] = (uint32_t)((uint16_t)(fOut0r - scratch0r)) |
                      (((uint16_t)(fOut0i - scratch0i)) << 16);

            // C_ADDTO(*Fout, t);
            Fout[0] = (uint32_t)((uint16_t)(fOut0r + scratch0r)) |
                      (((uint16_t)(fOut0i + scratch0i)) << 16);

            ++Fout;
        }
    }
}

void eff_fft4(eff_fft_cpx* src, eff_fft_cpx* dst) {
    // Computing log2(size)
    uint32_t log2Size = -1;
    uint32_t sizeShift = FFT_SIZE;
    while (sizeShift) {
        sizeShift >>= 1;
        log2Size++;
    }

    int shiftAmt = 32 - log2Size;

    // Copying and scrambling input
    eff_fft_init_dst(dst, src, FFT_SIZE, shiftAmt);

    // Running butterflies
    int i = 0;
    int maxStride = FFT_SIZE;
    for (int m = 1; m < maxStride; m *= 2) {
        int scheduleLen = maxStride / (m * 2);
        eff_fft_run_layer(dst, twiddleSchedule[i], m * 2, scheduleLen);
        i++;
    }
}

// A dummy function that triggers a bitstream miss and tile power up
// for each fabric function so that the main function doesn't need to
// be run twice
void eff_fft4_warmup() {
    int32_t fftSize = 0;
    stop_propagation_i32(&fftSize, 1);

    eff_fft_init_dst(0, 0, fftSize, 0);
    eff_fft_run_layer(0, 0, 0, 0);
}

#endif
