# MLPerfTiny Anomaly Detection

`mlperftiny_ad` is the MLPerfTiny anomaly-detection example for E1/E1x. It runs
the quantized `ad01_int8.tflite` autoencoder on a deterministic 640-sample int8
input vector and validates the reconstructed 640-sample output.

## What This App Demonstrates

Anomaly detection is a small neural network workload with a relatively compact
input and output tensor. It is useful for validating the end-to-end EFF model
flow because the output is a full reconstruction tensor rather than only a small
classification vector.

This app demonstrates:

- importing a quantized TFLite model with `eff-import`;
- generating separate fabric and scalar MLIR outputs;
- running the model inside an EFF profiling region;
- comparing every output element against a fixed expected tensor.

## Model and Tensor Contract

| Item | Value |
|------|-------|
| Model file | `ad01_int8.tflite` |
| Model source | MLCommons Tiny anomaly detection |
| Input tensor | `int8_t input[1][640]` |
| Output tensor | `int8_t output[1][640]` |
| Expected output | 640-element int8 reconstruction in `eff/benchmark.c` |
| Build targets | `fabric`, `scalar` |
| Architectures | `e1x`, `e1` |
| Operation count tag | `529416` |

## Build Flow

`CMakeLists.txt` locates `eff-import` next to the configured EFF C compiler and
uses it to generate:

- `ad.tosa.fabric.mlir` for the fabric target;
- `ad.tosa.scalar.mlir` for the scalar target.

The application sources are:

- `main.c`, which calls `benchmarkModel()`;
- `eff/benchmark.c`, which owns the input data, expected output, validation, and
  profiling loop;
- generated MLIR from the selected TFLite model.

## Running and Validation

Build from the parent `app_examples` folder:

```sh
mkdir bld
cd bld
cmake -G Ninja ..
ninja mlperftiny_ad_fabric
```

The app enters the `model` profile region, runs `run_model()` for the configured
profiling iteration count, and then compares all 640 output values. A passing
run prints:

```text
[mlperftiny_ad] PASS
```

Any output mismatch prints `FAIL` and returns a non-zero exit code.

## Customizing the Example

To try a different anomaly-detection model, replace `ad01_int8.tflite` and keep
the same input and output tensor contract, or update `eff/benchmark.c` to match
the new model shape. Regenerate the expected output from a known-good reference
run before using the app for validation.
