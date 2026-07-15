# Dense Matrix-Vector Multiply, Vectorized (DMV-vec)

A vectorized version of the dense matrix-vector multiply (`z = A · b`) for the Electron E1 general-purpose processor. It computes the same result as the base `dmv` example but packs the arithmetic into SIMD dot products.

## What's Different

This variant keeps the algorithm of the base kernel and changes how the inner loop is computed:

- It packs two 16-bit values per lane and uses the Fabric's SIMD dot-product operation to compute two multiply-accumulates at once, stepping the inner loop by 2.
- It processes 4 output rows per outer iteration, so four dot products advance together.
- Inputs are cast from 32-bit to 16-bit inline, which assumes the values fit in 16 bits.
- It assumes the row count is divisible by 4 and the column count is divisible by 2.

## Why EFF Hardware Performs Well

The Electron E1 runs this kernel on the Fabric architecture, a spatial dataflow design, so the effcc Compiler maps it onto the Fabric as a dataflow graph where operations fire as their inputs arrive. The SIMD dot-product operation doubles the arithmetic done per step, four rows accumulate in parallel, and partial sums stay resident on the Fabric while loads overlap with computation. The result is higher throughput per unit of energy than a scalar version of the same loop.

## Configurable Parameters

| Definition | Default | Effect |
|---|---|---|
| `NUM_ITERATIONS` | `1` | How many times the kernel runs. Increase it to average out noise when benchmarking. |
| `MAT_REF_SIZE` | `32` | The matrix dimension (both rows and columns), from the shared `mat.h`. It sets the problem size and must stay divisible by 4 for the row unrolling to be valid. |

> Full explanation of the algorithm, math, and applications:
> https://github.com/EfficientComputer/e1x_examples/tree/main/app_examples/dmv
