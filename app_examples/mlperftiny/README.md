# MLPerfTiny Benchmarks

This folder contains E1/E1x examples for the MLPerfTiny benchmark suite. Each
app builds a quantized TensorFlow Lite model with the EFF SDK, converts the
model through `eff-import`, runs a deterministic test input, and prints a
`[mlperftiny_<name>] PASS` or `FAIL` result.

Source snapshot: `EfficientComputer/apps` `main` at commit `46301721a6da`.

Upstream benchmark source: [mlcommons/tiny](https://github.com/mlcommons/tiny).
Pre-trained models were downloaded on April 6, 2026.

## Included Apps

| App | MLPerfTiny task | Active model | Input | Output | Build targets |
|-----|-----------------|--------------|-------|--------|---------------|
| `mlperftiny_ad` | Anomaly detection | `ad01_int8.tflite` | `[1][640]` int8 | `[1][640]` int8 reconstruction | `fabric`, `scalar` |
| `mlperftiny_ic` | Image classification | `models/v1.1/pretrainedResnet_quant.tflite` | `[1][32][32][3]` int8 image | `[1][10]` int8 logits | `fabric`, `optimized_fabric`, `scalar` |
| `mlperftiny_kws` | Keyword spotting | `kws_ref_model.tflite` | `[1][49][10][1]` int8 features | `[1][12]` int8 logits | `fabric`, `optimized_fabric`, `scalar` |
| `mlperftiny_sww` | Streaming wake word | `str_ww_ref_model.tflite` | 16000-sample int16 audio, then `[1][30][1][40]` int8 features | `[1][3]` int8 logits | `fabric`, `optimized_fabric`, `scalar` |
| `mlperftiny_vww` | Visual wake words | `vww_96_int8.tflite` | `[1][96][96][3]` int8 image | `[1][2]` int8 logits | `fabric`, `optimized_fabric`, `scalar` |

## Model Provenance

| Benchmark | Model file | Upstream commit | MD5 |
|-----------|------------|-----------------|-----|
| AD | [ad01_int8.tflite](https://github.com/mlcommons/tiny/blob/master/benchmark/training/anomaly_detection/trained_models/ad01_int8.tflite) | `bceb91c` | `361fa1b1b871e2068b2ab38d9805ef56` |
| IC v1.1 | [pretrainedResnet_quant.tflite](https://github.com/mlcommons/tiny/blob/360bf095/benchmark/training/image_classification/trained_models/pretrainedResnet_quant.tflite) | `360bf09` | `2d6dd48722471313e4c4528249205ae3` |
| IC v1.4-pre | [pretrainedResnet_quant.tflite](https://github.com/mlcommons/tiny/blob/eb78d0e/benchmark/training/image_classification/trained_models/pretrainedResnet_quant.tflite) | `eb78d0e` | `bd6cec63e4337ac66fb5ed5594b0df48` |
| KWS | [kws_ref_model.tflite](https://github.com/mlcommons/tiny/blob/master/benchmark/training/keyword_spotting/trained_models/kws_ref_model.tflite) | `bceb91c` | `618aeb155673dae92fa6a7f26608add5` |
| SWW | [str_ww_ref_model.tflite](https://github.com/mlcommons/tiny/blob/master/benchmark/training/streaming_wakeword/trained_models/str_ww_ref_model.tflite) | `904de6f` | `9643512e0b8a0b6e643e70c2fec496ee` |
| VWW | [vww_96_int8.tflite](https://github.com/mlcommons/tiny/blob/master/benchmark/training/visual_wake_words/trained_models/vww_96_int8.tflite) | `bceb91c` | `f0b011416abee0343a5d130cb1f4c18f` |

The image-classification folder includes both the v1.1 and v1.4-pre ResNet
models. The CMake file currently selects the v1.1 model; the v1.4-pre path is
kept next to it so developers can switch model revisions intentionally.

## Building

Build from `app_examples` using the same flow as the other E1x examples:

```sh
mkdir bld
cd bld
cmake -G Ninja .. -DEFF_STDIO_PORT=3
ninja
```

The MLPerfTiny apps are added through `app_examples/mlperftiny/CMakeLists.txt`,
which includes the five app subdirectories. Each app CMake file invokes
`eff-import` from the selected EFF compiler toolchain to translate the active
`.tflite` model into generated TOSA MLIR for the requested target.

## Running and Validation

Each app embeds a deterministic input tensor or sample, runs the model for the
configured number of profiling iterations, and compares the final output tensor
against an expected int8 result. A successful run prints:

```text
[mlperftiny_<name>] PASS
```

The examples are intended for normal EFF SDK simulator, scalar, and fabric
workflows. Optimized variants use custom rule files and hand-written kernels
where present.

## Folder Structure

Each app follows the same layout:

| Path | Purpose |
|------|---------|
| `CMakeLists.txt` | Selects the model, generated MLIR, targets, and sources. |
| `main.c` | Minimal entry point that calls `benchmarkModel()`. |
| `eff/benchmark.c` | Deterministic input, output validation, and profiling loop. |
| `eff/optimized_kernels*.c` | E1/E1x optimized kernels for apps with optimized builds. |
| `*.mlir`, `*.pdl` | Hand-written MLIR helpers and custom TOSA rewrite rules. |
| `arm/` | Reference ARM/CMSIS or Ethos-U ports retained from the source tree; not used by the E1x CMake targets. |

## TOSA Conversion

The app builds call `eff-import` automatically. For manual inspection, a
standalone TFLite to TOSA conversion can be attempted with:

```sh
pip install tosa-converter-for-tflite
tosa-converter-for-tflite --text <input.tflite> -o <output.tosa>
```

Manual output should be compared with `eff-import` output before it is used as a
replacement in these examples.
