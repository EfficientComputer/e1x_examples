// clang-format off <<AUTOBENCH>> efficient_e1 efficient_e1x ambiq_apollo4p nxp_lpc55s69 renesas_ra8m1 ambiq_apollo510 renesas_ra8p1 clang-format on

#include <limits.h>
#include <stdio.h>
#include <string.h>

#include "fft4k.h"
#include "mel_filterbank.h"
#include "hamming_coeffs.h"
#include "featurize.h"

// Implements a pre-emphasis (high-pass) filter
// Hard-coded single FIR tap of 31/32
__efficient__ void preemphasis(const int16_t *input, int16_t *output) {
    output[0] = input[0];
    for (size_t i = 1; i < WAV_NUM_SAMPLES; i++) {
        int32_t sample = (int32_t)input[i] - ((int32_t)input[i - 1] * 31 / 32);

        // Clip
        if (sample > INT16_MAX) {
            sample = INT16_MAX;
        } else if (sample < INT16_MIN) {
            sample = INT16_MIN;
        }

        output[i] = sample;
    }
}

__efficient__ void hamming_window(const int16_t *signal, const int16_t *coeffs,
                                  int16_t *out) {
    for (int f = 0; f < MEL_SPECTROGRAM_NUM_FRAMES; f++) {
        int inFrameStart = f * STFT_FRAME_STRIDE;
        int outFrameStart = f * STFT_FRAME_LENGTH;
        for (int i = 0; i < STFT_FRAME_LENGTH; i++) {
            int32_t y =
                (int32_t)signal[inFrameStart + i] * (int32_t)coeffs[i] / 32768;
            out[outFrameStart + i] = (int16_t)y;
        }
    }
}

__efficient__ void prepare_fft_inputs(
    eff_fft_cpx fft_arr[MEL_SPECTROGRAM_NUM_FRAMES][STFT_FRAME_LENGTH],
    int16_t real_input[STFT_FRAME_LENGTH * MEL_SPECTROGRAM_NUM_FRAMES]) {
    for (int frame = 0; frame < MEL_SPECTROGRAM_NUM_FRAMES; frame++) {
        int frameStart = frame * STFT_FRAME_LENGTH;
        for (int sample = 0; sample < STFT_FRAME_LENGTH; sample++) {
            fft_arr[frame][sample].r = real_input[frameStart + sample];
            fft_arr[frame][sample].i = 0;
        }
    }
}

// Mel filterbank is a bunch of triangular filters,
// basically a subproblem of dense-sparse matmul
// TODO: This currently fails on fabric because of int16_t -> int64_t promotion
// bug
__efficient__ void melspec_warp(
    int32_t power_spectrum[MEL_SPECTROGRAM_NUM_FRAMES][POWSPEC_LENGTH],
    int64_t mel_spectrum[MEL_SPECTROGRAM_NUM_FRAMES][MEL_SPECTROGRAM_NUM_BINS],
    const MelFilter *mel_filter_metadata, const int16_t *mel_weights) {
    for (int frame = 0; frame < MEL_SPECTROGRAM_NUM_FRAMES; frame++)  // 30
    {
        for (int bin = 0; bin < MEL_SPECTROGRAM_NUM_BINS; bin++)  // 40
        {
            MelFilter filter = mel_filter_metadata[bin];
            for (int i = filter.location; i < filter.location + filter.length;
                 i++) {
                mel_spectrum[frame][bin] +=
                    (int64_t)power_spectrum[frame][i] *
                    mel_weights[filter.offset + (i - filter.location)];
            }
        }
    }
}

__efficient__ void get_power_spectra(eff_fft_cpx fft[][STFT_FRAME_LENGTH],
                                     int32_t spectra[][POWSPEC_LENGTH]) {
    for (int frame = 0; frame < MEL_SPECTROGRAM_NUM_FRAMES; frame++) {
        for (int i = 0; i < POWSPEC_LENGTH; i++) {
            int16_t real = (int16_t)fft[frame][i].r;
            int16_t imag = (int16_t)fft[frame][i].i;
            // No risk of overflow due to FFT normalization
            int32_t sum = (real * real) + (imag * imag);
            // magnitude is sqrt(sum), but it gets squared again for powspec
            spectra[frame][i] = sum / STFT_FRAME_LENGTH;
        }
    }
}

// TODO: a lot of opportunities for early exit in this function but will have to
// explore
__efficient__ void log_normalize(
    int64_t input[MEL_SPECTROGRAM_NUM_FRAMES][MEL_SPECTROGRAM_NUM_BINS],
    int8_t *output) {
    for (int frame = 0; frame < MEL_SPECTROGRAM_NUM_FRAMES; frame++) {
        for (int bin = 0; bin < MEL_SPECTROGRAM_NUM_BINS; bin++) {
            int64_t m = input[frame][bin];
            int32_t log10_q8, log2_q8, frac, norm, quantized;
            if (m < 1) {
                log10_q8 = 0;
            } else {
                // Logarithm function manually inlined
                int n = 63 - __builtin_clzll((uint64_t)m);
                if (n >= 8) {
                    frac = (uint32_t)((uint64_t)m >> (n - 8)) & 0xFF;
                } else {
                    frac = (uint32_t)((uint64_t)m << (8 - n)) & 0xFF;
                }
                log2_q8 = (n << 8) | frac;
                log10_q8 = (log2_q8 * 771) >> 8;
            }

            norm = (log10_q8 - 5955) / 64;
            if (norm < 0) {
                norm = 0;
            } else if (norm > 256) {
                norm = 256;
            }

            quantized = ((norm * MODEL_SCALE_FACTOR_RECPR_INT + 128) >> 8) +
                        MODEL_ZERO_POINT;
            if (quantized < INT8_MIN) {
                quantized = INT8_MIN;
            } else if (quantized > INT8_MAX) {
                quantized = INT8_MAX;
            }
            output[frame * MEL_SPECTROGRAM_NUM_BINS + bin] = (int8_t)quantized;
        }
    }
}

int16_t filtered[WAV_NUM_SAMPLES];
int16_t afterHamming[STFT_FRAME_LENGTH * MEL_SPECTROGRAM_NUM_FRAMES];
int32_t powspec[MEL_SPECTROGRAM_NUM_FRAMES][POWSPEC_LENGTH];
int64_t melspec[MEL_SPECTROGRAM_NUM_FRAMES][MEL_SPECTROGRAM_NUM_BINS];

// Apollo4p can't fit FFT buffers in TCM. Put them in shared SRAM
#ifdef AM_PART_APOLLO4P
__attribute__((section(".shared"))) eff_fft_cpx
    fft_in[MEL_SPECTROGRAM_NUM_FRAMES][STFT_FRAME_LENGTH];
__attribute__((section(".shared"))) eff_fft_cpx
    fft_out[MEL_SPECTROGRAM_NUM_FRAMES][STFT_FRAME_LENGTH];
#else
eff_fft_cpx fft_in[MEL_SPECTROGRAM_NUM_FRAMES][STFT_FRAME_LENGTH];
eff_fft_cpx fft_out[MEL_SPECTROGRAM_NUM_FRAMES][STFT_FRAME_LENGTH];
#endif

void featurize(const int16_t *rawData, int8_t *featurized) {
    // Q1.15 -> holds [-1, 1) with 15 bit precision
    preemphasis(rawData, filtered);

    // Short time fourier transform
    // All intermediate results are buffered so each step of
    // the STFT can be batched

    hamming_window(filtered, hamming_coeffs, afterHamming);
    prepare_fft_inputs(fft_in, afterHamming);

    // NOP for Efficient devices, FFT init for others
    fft_init();

    for (int frame = 0; frame < MEL_SPECTROGRAM_NUM_FRAMES; frame++) {
        eff_fft4(fft_in[frame], fft_out[frame]);
    }

    get_power_spectra(fft_out, powspec);

    // Must zero out melspec before accumulation
    memset(melspec, 0, sizeof(melspec));
    melspec_warp(powspec, melspec, mel_filters, mel_filter_weights);

    log_normalize(melspec, featurized);
}
