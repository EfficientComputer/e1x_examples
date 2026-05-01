# MLPerfTiny Image Classification

`mlperftiny_ic` is the MLPerfTiny image-classification example for E1/E1x. It
runs a quantized ResNet model on a deterministic 32 x 32 RGB input image and
validates the 10-class output tensor.

## What This App Demonstrates

Image classification is the largest MLPerfTiny workload in this suite by the
operation-count tag. It exercises convolution-heavy inference, custom TOSA
rewrite rules, and E1/E1x optimized kernels.

This app demonstrates:

- selecting an MLPerfTiny image-classification TFLite model;
- building scalar, fabric, and optimized fabric variants;
- adding hand-written MLIR helpers for common low-level operations;
- switching between generic `run_model()` and optimized `run_model_opt()`;
- validating logits against a deterministic expected output.

## Model and Tensor Contract

| Item | Value |
|------|-------|
| Active model file | `models/v1.1/pretrainedResnet_quant.tflite` |
| Included alternate model | `models/v1.4-pre/pretrainedResnet_quant.tflite` |
| Model source | MLCommons Tiny image classification |
| Input tensor | `int8_t input[1][32][32][3]` |
| Output tensor | `int8_t output[1][10]` |
| Default sample | `CAT` input enabled in `eff/benchmark.c` |
| Expected output | `{-128, -128, -128, 122, -126, -125, -128, -128, -128, -128}` |
| Build targets | `fabric`, `optimized_fabric`, `scalar` |
| Architectures | `e1x`, `e1` |
| Operation count tag | `25113821` |

## Build Flow

`CMakeLists.txt` uses `eff-import` to generate three model forms:

- `ic.tosa.fabric.mlir` for the default fabric path;
- `ic.tosa.opt.fabric.mlir` for the optimized fabric path with custom kernels;
- `ic.tosa.scalar.mlir` for the scalar path.

The optimized path passes:

```text
--use-custom-kernels --custom-rule-path=eff/tosa_optimized_rules_${EFF_ARCH}.pdl --entrypoint-name=run_model_opt
```

Hand-written MLIR helpers are included from:

- `sdot_vec4xi8.mlir`;
- `add_tree.mlir`;
- `load_guard.mlir`.

Architecture-specific optimized kernels live in:

- `eff/optimized_kernels_e1.c`;
- `eff/optimized_kernels_e1x.c`.

## Running and Validation

Build from the parent `app_examples` folder:

```sh
mkdir bld
cd bld
cmake -G Ninja ..
ninja mlperftiny_ic_fabric
ninja mlperftiny_ic_optimized_fabric
```

The generic target calls `run_model()`. The optimized target defines
`USE_PRAGMAS`, which selects the optimized generated entry point and companion
kernels. A passing run prints:

```text
[mlperftiny_ic] PASS
```

Any mismatch in the 10 output logits prints `FAIL` and returns a non-zero exit
code.

## Customizing the Example

To test the included v1.4-pre model, change `MODEL_INPUT_PATH` in
`CMakeLists.txt` from the v1.1 path to the commented v1.4-pre path, then update
the expected output in `eff/benchmark.c` from a reference run. To test a new
image, replace the embedded input tensor and regenerate the expected logits.
