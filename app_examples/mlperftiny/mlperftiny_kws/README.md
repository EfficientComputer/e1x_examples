# MLPerfTiny Keyword Spotting

`mlperftiny_kws` is the MLPerfTiny keyword-spotting example for E1/E1x. It runs
the quantized keyword-spotting model on precomputed audio features and validates
the 12-class output tensor.

## What This App Demonstrates

Keyword spotting is a compact audio ML workload that stresses depthwise and
pointwise convolution patterns. This implementation includes both generic and
optimized paths so developers can compare generated model code with custom
E1/E1x kernels.

This app demonstrates:

- importing `kws_ref_model.tflite` with `eff-import`;
- building fabric, optimized fabric, and scalar targets;
- using custom TOSA rewrite rules for optimized inference;
- selecting E1 or E1x optimized kernels based on `EFF_ARCH`;
- validating a 12-class output vector from deterministic input features.

## Model and Tensor Contract

| Item | Value |
|------|-------|
| Model file | `kws_ref_model.tflite` |
| Model source | MLCommons Tiny keyword spotting |
| Input tensor | `int8_t input[1][49][10][1]` |
| Output tensor | `int8_t output[1][12]` |
| Expected output | `{-128, -128, -128, -128, -128, -128, -128, 127, -128, -128, -128, -128}` |
| Build targets | `fabric`, `optimized_fabric`, `scalar` |
| Architectures | `e1x`, `e1` |
| Operation count tag | `5393635` |

## Build Flow

`CMakeLists.txt` uses `eff-import` to generate:

- `kws.tosa.fabric.mlir` for fabric;
- `kws.tosa.opt.fabric.mlir` for optimized fabric with custom kernels;
- `kws.tosa.scalar.mlir` for scalar.

The optimized model generation uses:

```text
--use-custom-kernels --custom-rule-path=eff/tosa_optimized_rules_${EFF_ARCH}.pdl --entrypoint-name=run_model_opt
```

The build chooses the optimized kernel source by architecture:

- `eff/optimized_kernels_e1x.c` for E1x;
- `eff/optimized_kernels_e1.c` for E1.

The optimized target defines both `USE_PRAGMAS` and
`EFF_BLD_HAND_OPTIMIZED=1`, which makes `eff/benchmark.c` call
`run_model_opt()`.

## Running and Validation

Build from the parent `app_examples` folder:

```sh
mkdir bld
cd bld
cmake -G Ninja ..
ninja mlperftiny_kws_fabric
ninja mlperftiny_kws_optimized_fabric
```

The app enters the `model` profile region, runs the selected model entry point,
and compares all 12 logits against the expected output. A passing run prints:

```text
[mlperftiny_kws] PASS
```

Any mismatch prints `FAIL` and returns a non-zero exit code.

## Customizing the Example

To run a new feature sample, replace `input_arr` in `eff/benchmark.c` with a
49 x 10 x 1 int8 feature tensor and regenerate the expected 12-class output.
When changing the TFLite model, keep the tensor shape contract aligned with the
model or update the C types and validation loop accordingly.
