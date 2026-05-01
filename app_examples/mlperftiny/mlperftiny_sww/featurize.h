// clang-format off <<AUTOBENCH>> efficient_e1 efficient_e0 efficient_e1x ambiq_apollo4p nxp_lpc55s69 renesas_ra8m1 ambiq_apollo510 renesas_ra8p1 clang-format on

#pragma once
#ifdef __cplusplus
extern "C" {
#endif
#include <stdint.h>

#include "mel_filterbank.h"

#define WAV_NUM_SAMPLES 16000
#define MEL_SPECTROGRAM_NUM_BINS 40
#define MEL_SPECTROGRAM_NUM_FRAMES 30

#define STFT_FRAME_LENGTH 1024
#define STFT_FRAME_STRIDE 512

#define MODEL_SCALE_FACTOR 0.003701042616739869f
#define MODEL_SCALE_FACTOR_RECPR_INT ((int)(1.0f / MODEL_SCALE_FACTOR))
#define MODEL_ZERO_POINT -128

// Only track power spectrum of unique components of STFT
#define POWSPEC_LENGTH ((STFT_FRAME_LENGTH / 2) + 1)

void preemphasis(const int16_t *input, int16_t *output);
void hamming_window(const int16_t *signal, const int16_t *coeffs, int16_t *out);
void prepare_fft_inputs(
    eff_fft_cpx fft_arr[MEL_SPECTROGRAM_NUM_FRAMES][STFT_FRAME_LENGTH],
    int16_t real_input[STFT_FRAME_LENGTH * MEL_SPECTROGRAM_NUM_FRAMES]);
void melspec_warp(
    int32_t power_spectrum[MEL_SPECTROGRAM_NUM_FRAMES][POWSPEC_LENGTH],
    int64_t mel_spectrum[MEL_SPECTROGRAM_NUM_FRAMES][MEL_SPECTROGRAM_NUM_BINS],
    const MelFilter *mel_filter_metadata, const int16_t *mel_weights);
void log_normalize(
    int64_t input[MEL_SPECTROGRAM_NUM_FRAMES][MEL_SPECTROGRAM_NUM_BINS],
    int8_t *output);
void get_power_spectra(eff_fft_cpx fft[][STFT_FRAME_LENGTH],
                       int32_t spectra[][POWSPEC_LENGTH]);

void featurize(const int16_t *rawData, int8_t *featurized);

#ifdef __cplusplus
}
#endif
