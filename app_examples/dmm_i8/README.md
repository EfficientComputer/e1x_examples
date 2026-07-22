# Dense Matrix Multiply, int8 (DMM-i8)

An int8 version of the dense matrix multiply (`Z = A * B`) for the Electron E1 general-purpose processor. It computes the same operation as the base `dmm` example but on 8-bit integer inputs, with a quantization step applied to the results.

## What's Different

This variant keeps the matrix-multiply algorithm and changes the data type and output handling:

- Inputs `A` and `B` are 8-bit integers rather than 32-bit integers, so four values pack into each 32-bit word and are unpacked inside the inner loop.
- Products accumulate into wider 32-bit sums so the running total does not overflow.
- After the multiply-accumulate, each result is requantized back to 8 bits using a fixed-point scale: it is multiplied by an integer multiplier, rounded, and right-shifted by a configurable amount. This is the standard integer-only quantization used in low-precision inference.
- It keeps both the wide unscaled result and the final 8-bit output, and checks correctness with a CRC over the 8-bit output.

## Why EFF Hardware Performs Well

The Electron E1 runs this kernel on the Fabric architecture, a spatial dataflow design, so the effcc Compiler maps it onto the Fabric as a dataflow graph where operations fire as their inputs arrive. Packing four 8-bit values per word lets more multiply-accumulates advance per step while the tile accumulators stay resident on the Fabric, and the requantization multiply and shift fold directly into the same dataflow graph so the low-precision output is produced with very little extra energy.

## Configurable Parameters

| Definition | Default | Effect |
|---|---|---|
| `NUM_ITERATIONS` | `1` | How many times the kernel runs. Increase it to average out noise when benchmarking. |
| `multiplier` | `1` | The integer multiplier applied during requantization of each result. Changing it rescales the 8-bit output. |
| `shift` | `0` | The right-shift amount applied during requantization. It sets the fixed-point scale together with `multiplier`. |
| `MAT_REF_SIZE` | `32` | The matrix dimension (used for all of `n`, `m`, and `o`), from the shared `mat.h`. It sets the problem size and must stay divisible by 4 for the packing and tiling to be valid. |

> Full explanation of the algorithm, math, and applications:
> https://github.com/EfficientComputer/e1x_examples/tree/main/app_examples/dmm
