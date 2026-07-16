# 5x5 2D Convolution, Symmetric (Conv5x5-sym)

A variant of the 5x5 2D convolution for the Electron E1x general-purpose processor. It computes the same operation as the base `conv5x5` example, sliding a 5x5 filter across an input image, but exploits filter symmetry to cut the number of multiplies per output pixel.

## What's Different

This variant assumes the 5x5 filter is symmetric, so pixels at mirror-image positions in the window share the same weight:

- Input pixels that map to the same coefficient are added together first, then multiplied by that shared weight once, instead of multiplying each pixel separately.
- Because a symmetric 5x5 filter has 9 unique weights rather than 25, this reduces the multiply count per output pixel from 25 toward 9 while producing the same result.
- The filter is symmetrized at startup so the built-in or randomized weights satisfy the symmetry assumption.
- Work is done in tiles of output rows, and the code assumes the tile size divides the image dimension evenly.

## Why EFF Hardware Performs Well

The Electron E1x runs this kernel on the Fabric architecture, a spatial dataflow design, so the effcc Compiler maps it onto the Fabric as a dataflow graph where operations fire as their inputs arrive. Adding the symmetric pixel pairs before multiplying removes multiplies from the graph, the shared weights stay resident on the Fabric while the window slides, and partial sums accumulate on the Fabric while loads overlap with computation. The result is higher throughput per unit of energy than a version that multiplies every pixel independently.

## Configurable Parameters

| Definition         | Default | Effect                                                                                                                                                                   |
| ------------------ | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `NUM_ITERATIONS`   | `1`     | How many times the kernel runs. Increase it to average out noise when benchmarking.                                                                                      |
| `N`                | `64`    | The side length of the input image (both width and height), which sets the problem size.                                                                                 |
| `N_PAD`            | `N + 4` | The padded side length of the input buffer. The filter is 5 wide, so the buffer is padded by 4. Keep it as `N + 4`.                                                      |
| `RANDOMIZE_FILTER` | `1`     | When nonzero, the filter weights are overwritten with pseudo-random values at startup instead of the built-in weights. The weights are symmetrized afterward either way. |
| `RANGE`            | `8`     | The range of the pseudo-random filter weights, centered around zero.                                                                                                     |
| `TILE_ROWS`        | `N`     | The number of output rows processed per tile. It must divide `N` evenly.                                                                                                 |

> Full explanation of the algorithm, math, and applications:
> https://github.com/EfficientComputer/e1x_examples/tree/main/app_examples/conv5x5
