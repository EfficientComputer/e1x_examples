#define DEFINE_TWIDDLES
#include "fft.h"

#define FIXED_ROUND(x) \
    (fft_scalar)(((x) + (1 << (FRACBITS - 1))) >> FRACBITS)

#define SAMP_MAX 32767
#define FRACBITS 15

__efficient__ 
void fft_init_dst(fft_cpx* dst, fft_cpx* src,
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

__efficient__ 
void kiss_fft_run_layer(fft_cpx* data, int twiddleStart,
                                   int idxStride, int scheduleLen) {
    for (int i = 0; i < scheduleLen; i++) {
        uint32_t* Fout = (uint32_t*)(data + (i * idxStride));
        uint32_t m = idxStride / 4;

        uint32_t *tw1, *tw2, *tw3;
        uint32_t m2 = m + m;
        uint32_t m3 = m2 + m;

        tw3 = tw2 = tw1 = (uint32_t*)(twiddles + twiddleStart * 2);

        for (int32_t k = m; k > 0; k--) {
            // Equivalent to *Fout / 4;
            uint32_t Fout0 = Fout[0];
            int16_t Fout0r = (int16_t)(Fout0 & 0xFFFF);
            int16_t Fout0i = (int16_t)(Fout0 >> 16);
            int32_t fOut0r = FIXED_ROUND((int32_t)(Fout0r * (SAMP_MAX / 4)));
            int32_t fOut0i = FIXED_ROUND((int32_t)(Fout0i * (SAMP_MAX / 4)));
            
            // Equivalent to Fout[m] / 4;
            uint32_t Fout1 = Fout[m];
            int16_t Fout1r = (int16_t)(Fout1 & 0xFFFF);
            int16_t Fout1i = (int16_t)(Fout1 >> 16);
            int32_t fOut1r = FIXED_ROUND((int32_t)(Fout1r * (SAMP_MAX / 4)));
            int32_t fOut1i = FIXED_ROUND((int32_t)(Fout1i * (SAMP_MAX / 4)));
            
            // Equivalent to Fout[m2] / 4;
            uint32_t Fout2 = Fout[m2];
            int16_t Fout2r = (int16_t)(Fout2 & 0xFFFF);
            int16_t Fout2i = (int16_t)(Fout2 >> 16);
            int32_t fOut2r = FIXED_ROUND((int32_t)(Fout2r * (SAMP_MAX / 4)));
            int32_t fOut2i = FIXED_ROUND((int32_t)(Fout2i * (SAMP_MAX / 4)));
            
            // Equivalent to Fout[m3] / 4;
            uint32_t Fout3 = Fout[m3];
            int16_t Fout3r = (int16_t)(Fout3 & 0xFFFF);
            int16_t Fout3i = (int16_t)(Fout3 >> 16);
            int32_t fOut3r = FIXED_ROUND((int32_t)(Fout3r * (SAMP_MAX / 4)));
            int32_t fOut3i = FIXED_ROUND((int32_t)(Fout3i * (SAMP_MAX / 4)));

            // Equivalent to scratch[0] = Fout[m] * *tw1;
            uint32_t twiddle1 = *tw1;
            int16_t twiddle1r = (int16_t)(twiddle1 & 0xFFFF);
            int16_t twiddle1i = (int16_t)(twiddle1 >> 16);
            int32_t scratch0r =
                FIXED_ROUND((int32_t)(fOut1r * twiddle1r - fOut1i * twiddle1i));
            int32_t scratch0i =
                FIXED_ROUND((int32_t)(fOut1r * twiddle1i + fOut1i * twiddle1r));

            // Equivalent to scratch[1] = Fout[m2] * *tw2;
            uint32_t twiddle2 = *tw2;
            int16_t twiddle2r = (int16_t)(twiddle2 & 0xFFFF);
            int16_t twiddle2i = (int16_t)(twiddle2 >> 16);
            int32_t scratch1r =
                FIXED_ROUND((int32_t)(fOut2r * twiddle2r - fOut2i * twiddle2i));
            int32_t scratch1i =
                FIXED_ROUND((int32_t)(fOut2r * twiddle2i + fOut2i * twiddle2r));

            // Equivalent to scratch[2] = Fout[m3] * *tw3;
            uint32_t twiddle3 = *tw3;
            int16_t twiddle3r = (int16_t)(twiddle3 & 0xFFFF);
            int16_t twiddle3i = (int16_t)(twiddle3 >> 16);
            int32_t scratch2r =
                FIXED_ROUND((int32_t)(fOut3r * twiddle3r - fOut3i * twiddle3i));
            int32_t scratch2i =
                FIXED_ROUND((int32_t)(fOut3r * twiddle3i + fOut3i * twiddle3r));

            // Equivalent to scratch[5] = *Fout - scratch[1;
            int32_t scratch5r = fOut0r - scratch1r;
            int32_t scratch5i = fOut0i - scratch1i;

            // Equivalent to *Fout *  scratch[1];
            fOut0r = fOut0r + scratch1r;
            fOut0i = fOut0i + scratch1i;

            // Equivalent to scratch[3] = scratch[0] + scratch[2];
            int32_t scratch3r = scratch0r + scratch2r;
            int32_t scratch3i = scratch0i + scratch2i;

            // Equivalent to scratch[4] = scratch[0] - scratch[2];
            int32_t scratch4r = scratch0r - scratch2r;
            int32_t scratch4i = scratch0i - scratch2i;

            // Equivalent to Fout[m2] = *Fout - scratch[3];
            Fout[m2] = (uint32_t)((uint16_t)(fOut0r - scratch3r)) |
                       (((uint16_t)(fOut0i - scratch3i)) << 16);

            tw1 += 1;
            tw2 += 2;
            tw3 += 3;
            // Equivalent to *Fout + scratch[3];
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

void fft4(fft_cpx* src, fft_cpx* dst) {
    fft_init_dst(dst, src, FFT_SIZE);

    int i = 0;
    int maxStride = FFT_SIZE;
    for (int m = 1; m < maxStride; m *= 4) {
        int scheduleLen = maxStride / (m * 4);
        kiss_fft_run_layer(dst, twiddleSchedule[i], m * 4, scheduleLen);
        i++;
    }
}
