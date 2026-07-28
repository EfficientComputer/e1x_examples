# Dense Matrix-Vector Multiply, Floating-Point (DMV-fp)

This example is a floating-point version of the dense matrix-vector multiply (`z = A · b`) for the Electron E1x general-purpose processor. It computes the same operation as the base `dmv` example but on floating-point data instead of integers.

## What's Different

This variant keeps the matrix-vector algorithm and changes the data type and loop shape:

- It operates on floating-point values rather than 32-bit integers.
- The optimized build unrolls 8 output rows at a time and reuses each vector element across those rows, which cuts the number of loads in the inner loop.
- It uses a separate row stride (`MSTRIDE`) so matrix rows can be padded to a convenient width, independent of the logical column count.
- It includes both a correctness test (against a hand-computed reference) and a benchmark entry point.

## Why Efficient Hardware Performs Well

The Electron E1x runs this kernel on the Fabric architecture, a spatial dataflow design, so the effcc Compiler maps it onto the Fabric as a dataflow graph where operations fire as their inputs arrive. Eight row accumulators stay resident on the Fabric and advance in parallel, reusing each loaded vector element across all eight rows so memory traffic stays low while the floating-point multiply-accumulates overlap with the remaining loads.

## Configurable Parameters

| Definition       | Default | Effect                                                                                                          |
| ---------------- | ------- | --------------------------------------------------------------------------------------------------------------- |
| `NUM_ITERATIONS` | `1`     | This is how many times the benchmark kernel runs. Increase it to average out noise when benchmarking.                   |
| `M`              | `32`    | The number of matrix rows (output length)                                                                      |
| `N`              | `32`    | The vector length and number of matrix columns                                                                 |
| `MSTRIDE`        | `34`    | This is the row stride in elements. It must be at least `N`; the gap between `N` and `MSTRIDE` is padding between rows. |

> Full explanation of the algorithm, math, and applications:
> https://github.com/EfficientComputer/e1x_examples/tree/main/app_examples/dmv
