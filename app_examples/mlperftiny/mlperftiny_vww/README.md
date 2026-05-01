# MLPerfTiny Visual Wake Words

`mlperftiny_vww` is the MLPerfTiny visual-wake-words example for E1/E1x. It
runs the quantized `vww_96_int8.tflite` model on a deterministic 96 x 96 RGB
input image and validates the two-class output tensor.

## What This App Demonstrates

Visual wake words is an image-based binary classification workload. It is useful
for exercising convolutional inference with a larger image tensor than the
32 x 32 image-classification example while keeping the output contract simple.

This app demonstrates:

- importing the MLPerfTiny VWW TFLite model with `eff-import`;
- building scalar, fabric, and optimized fabric variants;
- using custom rewrite rules and hand-written MLIR helpers;
- validating a compact two-logit output tensor.

## Model and Tensor Contract

| Item | Value |
|------|-------|
| Model file | `vww_96_int8.tflite` |
| Model source | MLCommons Tiny visual wake words |
| Input tensor | `int8_t input[1][96][96][3]` |
| Output tensor | `int8_t output[1][2]` |
| Expected output | `{-114, 114}` |
| Build targets | `fabric`, `optimized_fabric`, `scalar` |
| Architectures | `e1x`, `e1` |
| Operation count tag | `15213445` |

## Build Flow

`CMakeLists.txt` uses `eff-import` to generate:

- `vww.tosa.fabric.mlir` for fabric;
- `vww.tosa.opt.fabric.mlir` for optimized fabric with custom kernels;
- `vww.tosa.scalar.mlir` for scalar.

The optimized path currently uses `eff/tosa_optimized_rules.mlir`.
Architecture-specific rule files can be introduced when E1 support is added for
that path.

Hand-written MLIR helpers are included from:

- `sdot_vec4xi8.mlir`;
- `add_tree.mlir`;
- `load_guard.mlir`.

The optimized kernel implementation is `eff/optimized_kernels.c`.

## Running and Validation

Build from the parent `app_examples` folder:

```sh
mkdir bld
cd bld
cmake -G Ninja ..
ninja mlperftiny_vww_fabric
ninja mlperftiny_vww_optimized_fabric
```

The app enters the `model` profile region, runs the selected model entry point,
and compares both output logits against the expected result. A passing run
prints:

```text
[mlperftiny_vww] PASS
```

Any mismatch prints `FAIL` and returns a non-zero exit code.

## Customizing the Example

To test a different visual input, replace `input_arr` in `eff/benchmark.c` with
a 96 x 96 x 3 int8 image tensor and regenerate the expected two-logit output.
When changing the model, update the TFLite path, any custom rewrite rules, and
the validation tensor together.
