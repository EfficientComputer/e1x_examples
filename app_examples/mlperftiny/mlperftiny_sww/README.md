# MLPerfTiny Streaming Wake Word

`mlperftiny_sww` is the MLPerfTiny streaming wake-word example for E1/E1x. It
runs a complete audio pipeline: a deterministic 16000-sample `marvin` waveform
is featurized into a mel-spectrogram tensor, then passed through the quantized
`str_ww_ref_model.tflite` classifier.

## What This App Demonstrates

Streaming wake word combines signal processing and neural network inference. It
is useful as an end-to-end audio example because the workload includes
preemphasis, framing, a Hamming window, FFTs, mel filtering, log normalization,
and model inference.

This app demonstrates:

- using checked-in generated DSP tables for deterministic builds;
- importing a quantized streaming wake-word TFLite model with `eff-import`;
- building scalar, fabric, and optimized fabric variants;
- using custom TOSA rewrite rules and E1/E1x optimized kernels;
- validating a three-class classifier output from a known audio sample.

## Model and Tensor Contract

| Item | Value |
|------|-------|
| Model file | `str_ww_ref_model.tflite` |
| Model source | MLCommons Tiny streaming wake word |
| Raw input | `int16_t input_data_marvin[16000]` |
| Feature tensor | `int8_t input[1][30][1][40]` |
| Output tensor | `int8_t output[1][3]` |
| Expected output | `{127, -128, -128}` |
| Build targets | `fabric`, `optimized_fabric`, `scalar` |
| Architectures | `e1x`, `e1` |
| Operation count tag | `1703877` |

## Feature Extraction

Key constants:

| Constant | Value |
|----------|-------|
| `WAV_NUM_SAMPLES` | `16000` |
| `STFT_FRAME_LENGTH` | `1024` |
| `STFT_FRAME_STRIDE` | `512` |
| `MEL_SPECTROGRAM_NUM_FRAMES` | `30` |
| `MEL_SPECTROGRAM_NUM_BINS` | `40` |
| `POWSPEC_LENGTH` | `513` |

The featurizer performs:

1. preemphasis over the raw audio;
2. frame extraction and Hamming windowing;
3. a 1024-point FFT per frame;
4. power spectrum calculation;
5. mel filterbank projection;
6. log normalization and int8 quantization.

The checked-in generated tables are:

- `hamming_coeffs.h`;
- `fft4_1024.h`;
- `mel_filterbank.h`.

The generator scripts are retained next to the generated headers so the tables
can be regenerated when the model or feature settings change.

## Approximate Operation Counts

| Stage | Approximate multiplies |
|-------|------------------------|
| `featurize()` | `877509` |
| model inference | `826368` |
| combined | `1703877` |

The FFT dominates the feature-extraction cost.

## Build Flow

`CMakeLists.txt` uses `eff-import` to generate:

- `sww.tosa.fabric.mlir` for fabric;
- `sww.tosa.opt.fabric.mlir` for optimized fabric with custom kernels;
- `sww.tosa.scalar.mlir` for scalar.

The optimized path uses `eff/tosa_optimized_rules_${EFF_ARCH}.pdl` and the
architecture-specific kernel file `eff/optimized_kernels_${EFF_ARCH}.c`.

## Running and Validation

Build from the parent `app_examples` folder:

```sh
mkdir bld
cd bld
cmake -G Ninja ..
ninja mlperftiny_sww_fabric
ninja mlperftiny_sww_optimized_fabric
```

The app enters the `model` profile region, runs featurization and inference for
the configured profiling iteration count, and compares all three output logits.
A passing run prints:

```text
[mlperftiny_sww] PASS
```

Any mismatch prints `FAIL` and returns a non-zero exit code.

## Customizing the Example

To test a different wake-word sample, replace the data in `marvin.h` with a
16000-sample signed 16-bit waveform and regenerate the expected three-class
output. If the feature settings change, regenerate the checked-in DSP tables
with the included Python scripts before rebuilding.
