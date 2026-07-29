# Finite Impulse Response Filter, Vectorized (FIR-vec)

This example is a vectorized version of the finite impulse response (FIR) filter for the Electron E1x general-purpose processor. It computes the same result as the base `fir` example but packs the arithmetic into 16-bit SIMD dot products.

## What's Different

This variant keeps the sliding-window filter algorithm and changes how each output sample is computed:

- The 16 filter taps are packed into 8 pairs of 16-bit values once per call and reused across every output sample.
- Each output sample is computed with 8 SIMD dot-product operations instead of 16 scalar multiply-accumulates, with each dot product handling two taps at once.
- The input window is held in a shift register of scalar values that advances by one sample per output, so overlapping windows reuse loaded data.
- It assumes 16 taps and that the values fit in 16 bits, and it checks correctness by comparing all outputs against an inline scalar reference.

## Why Efficient Hardware Performs Well

The Electron E1x runs this kernel on the Fabric architecture, a spatial dataflow design, so the effcc Compiler maps it onto the Fabric as a dataflow graph where operations fire as their inputs arrive. The SIMD dot-product operation halves the number of steps per output sample by processing two taps per lane, the packed taps stay resident on the Fabric and are reused for every sample, and the shift register keeps the overlapping input window on the Fabric so memory traffic stays low. The result is higher throughput per unit of energy than a scalar version of the same loop.

## Configurable Parameters

| Definition       | Default | Effect                                                                                                                                                                           |
| ---------------- | ------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `N`              | `64`    | This is the number of output samples produced, defined in `fir_vec.h`. A larger value means a longer output signal and more work per run.                                                |
| `W`              | `16`    | This is the number of filter taps, defined in `fir_vec.h`. The vectorized kernel assumes 16 taps packed into 8 pairs, so changing it also requires reworking the packing in `fir_vec.c`. |
| `NUM_ITERATIONS` | `1`     | This is how many times the kernel runs. Increase it to average out noise when benchmarking.                                                                                              |

The input signal `x` and the taps `w` are generated inline in `main.c` (`x[i] = i` and `w[i] = i - 8`), and the output is checked against an inline scalar reference, so no separate expected constant needs updating when the inputs change.

> Full explanation of the algorithm, math, and applications:
> https://github.com/EfficientComputer/e1x_examples/tree/main/app_examples/fir
