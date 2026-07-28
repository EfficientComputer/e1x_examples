# 5x5 2D Convolution, Vectorized (Conv5x5-vec)

This example is a vectorized version of the 5x5 2D convolution for the Electron E1x general-purpose processor. It computes the same result as the base `conv5x5` example, sliding a 5x5 filter across an input image, but packs the per-pixel arithmetic into SIMD dot products.

## What's Different

This variant keeps the algorithm of the base kernel and changes how the 25 multiply-accumulates per output pixel are computed:

- It packs two 16-bit values per lane and uses the Fabric's SIMD dot-product operation, so the 25 products become 12 packed dot products plus 1 scalar multiply.
- Filter weights are pre-packed into pairs once at the start, then reused for every output pixel.
- Inputs and weights are cast from 32-bit to 16-bit inline, which assumes the values fit in 16 bits.

## Why Efficient Hardware Performs Well

The Electron E1x runs this kernel on the Fabric architecture, a spatial dataflow design, so the effcc Compiler maps it onto the Fabric as a dataflow graph where operations fire as their inputs arrive. The SIMD dot-product operation doubles the arithmetic done per step, the packed filter weights stay resident on the Fabric while the window slides, and partial sums accumulate on the Fabric while loads overlap with computation. The result is higher throughput per unit of energy than a scalar version of the same convolution.

## Configurable Parameters

| Definition                 | Default | Effect                                                                                                                 |
| -------------------------- | ------- | ---------------------------------------------------------------------------------------------------------------------- |
| `NUM_ITERATIONS`           | `1`     | This is how many times the kernel runs. Increase it to average out noise when benchmarking.                                    |
| `N`                        | `64`    | The output width, which sets the problem size along one dimension                                                     |
| `M`                        | `64`    | The output height, which sets the problem size along the other dimension                                              |
| `CONV5x5_N_PAD`            | `N + 4` | This is the padded side length of the input buffer. The filter is 5 wide, so the buffer is padded by 4. Keep it as `N + 4`.    |
| `CONV5x5_RANDOMIZE_FILTER` | `1`     | When nonzero, the filter weights are overwritten with pseudo-random values at startup instead of the built-in weights. |
| `CONV5x5_RANGE`            | `8`     | This is the range of the pseudo-random filter weights, centered around zero. Keeping it small helps the values fit in 16 bits. |

> Full explanation of the algorithm, math, and applications:
> https://github.com/EfficientComputer/e1x_examples/tree/main/app_examples/conv5x5
