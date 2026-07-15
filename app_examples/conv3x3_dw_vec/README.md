# Depthwise 3x3 Convolution, Vectorized (Conv3x3-DW-vec)

A vectorized version of the depthwise 3x3 convolution for the Electron E1 general-purpose processor. It computes the same result as the base `conv3x3_dw` example, applying one independent 3x3 filter per channel, but packs the arithmetic into SIMD dot products.

## What's Different

This variant keeps the per-channel sliding-window algorithm of the base kernel and changes how the nine multiply-accumulates per output are computed:

- For each channel it pre-packs the nine int8 filter taps into four pairs (plus one leftover tap) and uses the Fabric's SIMD dot-product operation, replacing nine scalar multiply-accumulates per output with four paired dot products and one scalar multiply.
- Packing is done once per channel, hoisted above the pixel loops, and reused across every output position in that channel.
- Values are widened from int8 to 16-bit before the dot products, and the left and right edge columns keep their scalar handling to apply the SAME zero padding.

## Why EFF Hardware Performs Well

The Electron E1 runs this kernel on the Fabric architecture, a spatial dataflow design, so the effcc Compiler maps it onto the Fabric as a dataflow graph where operations fire as their inputs arrive. The channels are independent and their work proceeds in parallel, the SIMD dot-product operation doubles the arithmetic done per step, per-channel filter taps stay resident on the Fabric, and loads for the sliding window overlap with computation. The result is higher throughput per unit of energy than a scalar version of the same loop.

## Configurable Parameters

| Definition | Default | Effect |
|---|---|---|
| `NUM_ITERATIONS` | `1` | How many times the kernel runs. Increase it to average out noise when benchmarking. |
| `M` | `48` | The number of input rows (image height). |
| `N` | `48` | The number of input columns (image width). |
| `NCHANNELS` | `8` | The number of channels, each with its own independent 3x3 filter. |
| `KERNEL_DIM` | `3` | The filter side length. This kernel is specialized for 3x3 filters. |
| `STRIDE_H` | `1` | The vertical stride between output positions. |
| `STRIDE_W` | `1` | The horizontal stride between output positions. |
| `PAD_TOP` | `1` | Rows of zero padding added at the top. |
| `PAD_LEFT` | `1` | Columns of zero padding added on the left. |
| `PAD_RIGHT` | `1` | Columns of zero padding added on the right. |
| `PAD_BOTTOM` | `1` | Rows of zero padding added at the bottom. |
| `OH` | `M` | The output height. With SAME padding and stride 1 it equals the input height. |
| `OW` | `N` | The output width. With SAME padding and stride 1 it equals the input width. |

> Full explanation of the algorithm, math, and applications:
> https://github.com/EfficientComputer/e1x_examples/tree/main/app_examples/conv3x3_dw
