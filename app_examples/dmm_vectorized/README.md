# Dense Matrix Multiply, Vectorized (DMM-vec)

A vectorized version of the dense matrix multiply (`Z = A * B`) for the Electron E1x general-purpose processor. It computes the same result as the base `dmm` example but packs the arithmetic into 16-bit SIMD dot products.

## What's Different

This variant keeps the matrix-multiply algorithm and changes how the inner loop is computed:

- Inputs are cast from 32-bit to 16-bit integers, which assumes the values fit in 16 bits, and two 16-bit values are packed per lane.
- It uses a SIMD dot-product operation to compute two multiply-accumulates at once, stepping the shared inner dimension by 2.
- It computes a 4 by 4 tile of output elements per outer iteration, so sixteen accumulators advance together.
- The right input matrix is transposed up front in 4 by 4 blocks so both operands stream with regular access, and the accumulators are treated as independent so their memory-ordering constraints are relaxed.
- It assumes the row and column counts are divisible by 4 and the shared dimension is divisible by 2, and it checks correctness against an element-wise reference matmul rather than a fold.

## Why EFF Hardware Performs Well

The Electron E1x runs this kernel on the Fabric architecture, a spatial dataflow design, so the effcc Compiler maps it onto the Fabric as a dataflow graph where operations fire as their inputs arrive. The SIMD dot-product operation doubles the arithmetic done per step by processing two 16-bit values per lane, sixteen tile accumulators stay resident on the Fabric and advance in parallel, and loads overlap with computation. The result is higher throughput per unit of energy than a scalar version of the same loop.

## Configurable Parameters

| Definition       | Default | Effect                                                                                                                                                               |
| ---------------- | ------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `NUM_ITERATIONS` | `1`     | How many times the kernel runs. Increase it to average out noise when benchmarking.                                                                                  |
| `MAT_REF_SIZE`   | `32`    | The matrix dimension (used for all of `n`, `m`, and `o`), from the shared `mat.h`. It sets the problem size and must stay divisible by 4 for the tiling to be valid. |

> Full explanation of the algorithm, math, and applications:
> https://github.com/EfficientComputer/e1x_examples/tree/main/app_examples/dmm
