# 3x3 Convolution, Vectorized (Conv3x3-vec)

A vectorized version of the single-channel 3x3 2D convolution for the Electron E1x general-purpose processor. It computes the same result as the base `conv3x3` example but packs the arithmetic into SIMD dot products.

## What's Different

This variant keeps the sliding-window algorithm of the base kernel and changes how the nine multiply-accumulates per output are computed:

- It pre-packs the nine filter weights into four pairs (plus one leftover weight) and uses the Fabric's SIMD dot-product operation, replacing nine scalar multiply-accumulates per output with four paired dot products and one scalar multiply.
- Input values for each 3x3 window are packed into pairs the same way before being fed to the dot products.
- Values are cast from 32-bit to 16-bit inline, which assumes the input and filter values fit in 16 bits.

## Why EFF Hardware Performs Well

The Electron E1x runs this kernel on the Fabric architecture, a spatial dataflow design, so the effcc Compiler maps it onto the Fabric as a dataflow graph where operations fire as their inputs arrive. The SIMD dot-product operation doubles the arithmetic done per step, filter weights stay resident on the Fabric and are reused across every output, and loads for the sliding window overlap with computation. The result is higher throughput per unit of energy than a scalar version of the same loop.

## Configurable Parameters

| Definition       | Default | Effect                                                                              |
| ---------------- | ------- | ----------------------------------------------------------------------------------- |
| `NUM_ITERATIONS` | `1`     | How many times the kernel runs. Increase it to average out noise when benchmarking. |
| `M`              | `64`    | The number of input rows (image height).                                            |
| `N`              | `64`    | The number of input columns (image width).                                          |
| `INSTRIDE`       | `N`     | The row stride of the input buffer in elements.                                     |
| `OUTSTRIDE`      | `N`     | The row stride of the output buffer in elements.                                    |

> Full explanation of the algorithm, math, and applications:
> https://github.com/EfficientComputer/e1x_examples/tree/main/app_examples/conv3x3
