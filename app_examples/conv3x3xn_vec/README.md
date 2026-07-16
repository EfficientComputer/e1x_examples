# Multi-Channel 3x3 Convolution, Vectorized (Conv3x3xN-vec)

A vectorized version of the multi-channel 3x3 convolution for the Electron E1x general-purpose processor. It computes the same result as the base `conv3x3xn` example, convolving a 3x3 filter over every input channel and summing across channels, but packs the arithmetic into SIMD dot products.

## What's Different

This variant keeps the algorithm of the base kernel and changes how the nine multiply-accumulates per channel are computed:

- For each channel it pre-packs the nine filter taps into four pairs (plus one leftover tap) and uses the Fabric's SIMD dot-product operation, replacing nine scalar multiply-accumulates per output with four paired dot products and one scalar multiply.
- The channel loop is outermost, so per-channel weight packing happens once per channel rather than once per output pixel.
- The output is zero-initialized and each channel accumulates into it with a running sum. The channel loop stays sequential so cross-channel writes to the same output location keep their ordering.
- Values are cast from 32-bit to 16-bit inline before the dot products, which assumes the values fit in 16 bits.

## Why EFF Hardware Performs Well

The Electron E1x runs this kernel on the Fabric architecture, a spatial dataflow design, so the effcc Compiler maps it onto the Fabric as a dataflow graph where operations fire as their inputs arrive. The SIMD dot-product operation doubles the arithmetic done per step, per-channel filter taps stay resident on the Fabric, cross-channel partial sums accumulate on the Fabric, and loads for the sliding window overlap with computation. The result is higher throughput per unit of energy than a scalar version of the same loop.

## Configurable Parameters

| Definition       | Default              | Effect                                                                                                                                 |
| ---------------- | -------------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| `NUM_ITERATIONS` | `1`                  | How many times the kernel runs. Increase it to average out noise when benchmarking.                                                    |
| `M`              | `16`                 | The number of input rows (image height). Keeping it a power of two helps performance.                                                  |
| `N`              | `16`                 | The number of input columns (image width). Keeping it a power of two helps performance.                                                |
| `KERNEL_DIM`     | `3`                  | The filter side length. This kernel is specialized for 3x3 filters.                                                                    |
| `NCHANNELS`      | `4`                  | The number of input channels summed together for each output value.                                                                    |
| `FASTPATH`       | `1`                  | When enabled, adds padding columns to the input and output buffers so the store does not need to be gated, which improves performance. |
| `INSTRIDE`       | `N + 1`              | The row stride of the input buffer in elements. With `FASTPATH` enabled it includes one padding column to avoid memory bank conflicts. |
| `OUTSTRIDE`      | `N + KERNEL_DIM - 1` | The row stride of the output buffer in elements. With `FASTPATH` enabled it includes the leading padding columns.                      |

> Full explanation of the algorithm, math, and applications:
> https://github.com/EfficientComputer/e1x_examples/tree/main/app_examples/conv3x3xn
