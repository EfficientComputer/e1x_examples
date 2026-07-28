# 5x5 2D Convolution, Symmetric and Vectorized (Conv5x5-sym-vec)

This example is a vectorized version of the symmetric 5x5 2D convolution for the Electron E1x general-purpose processor. It computes the same operation as the base `conv5x5sym` example, sliding a symmetric 5x5 filter across an input image, but packs the reduced set of multiplies into SIMD dot products.

## What's Different

This variant combines the symmetry trick of the base kernel with vectorized arithmetic:

- As in the symmetric base kernel, input pixels that share a filter weight are summed first so each of the 9 unique weights is multiplied only once.
- Those 9 unique multiplies are then packed into 4 SIMD dot products (two 16-bit values per lane) plus 1 scalar multiply, using the Fabric's SIMD dot-product operation.
- The 9 unique filter weights are pre-packed into pairs once at the start and reused for every output pixel.
- Because up to four input pixels are summed before multiplying, the inputs must fit in a 13-bit range so the intermediate sums stay within 16 bits.

## Why Efficient Hardware Performs Well

The Electron E1x runs this kernel on the Fabric architecture, a spatial dataflow design, so the effcc Compiler maps it onto the Fabric as a dataflow graph where operations fire as their inputs arrive. Summing the symmetric pixel pairs first removes multiplies, the SIMD dot-product operation then doubles the arithmetic done per remaining step, and the packed weights and partial sums stay resident on the Fabric while loads overlap with computation. The result is higher throughput per unit of energy than either the scalar symmetric version or a plain vectorized convolution.

## Configurable Parameters

| Definition         | Default | Effect                                                                                                                                                                   |
| ------------------ | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `NUM_ITERATIONS`   | `1`     | This is how many times the kernel runs. Increase it to average out noise when benchmarking.                                                                                      |
| `N`                | `64`    | The side length of the input image (both width and height), which sets the problem size                                                                                 |
| `N_PAD`            | `N + 4` | This is the padded side length of the input buffer. The filter is 5 wide, so the buffer is padded by 4. Keep it as `N + 4`.                                                      |
| `RANDOMIZE_FILTER` | `1`     | When nonzero, the filter weights are overwritten with pseudo-random values at startup instead of the built-in weights. The weights are symmetrized afterward either way. |
| `RANGE`            | `8`     | This is the range of the pseudo-random filter weights, centered around zero. Keeping it small helps the summed inputs stay within 16 bits.                                       |
| `TILE_ROWS`        | `N`     | This is the number of output rows processed per tile. It must divide `N` evenly.                                                                                                 |

> Full explanation of the algorithm, math, and applications:
> https://github.com/EfficientComputer/e1x_examples/tree/main/app_examples/conv5x5sym
