# Fast Fourier Transform (FFT)

This example computes a fixed-point Fast Fourier Transform on the Electron E1x general-purpose processor. It is a radix-4, Q15 FFT that transforms a block of time-domain samples into the frequency domain, a core building block for digital signal processing.

---

## Overview

### What is an FFT?

A Fast Fourier Transform is an efficient algorithm for computing the Discrete Fourier Transform (DFT), which converts a signal from the time domain into the frequency domain. This example uses a radix-4 decomposition, which processes four points at a time and reduces the work from O(N^2) for a direct DFT to O(N log N). The input and output are stored as complex numbers in Q15 fixed-point format, and the output bins are produced in bit-reversed order.

### Mathematical Definition

    X[k] = Σ x[n] * e^(-j*2*pi*k*n/N)     for n = 0 .. N-1

Where:

- `x[n]` is the input signal in the time domain (complex, Q15)
- `X[k]` is the output in the frequency domain (complex, Q15)
- `k` is the frequency bin index
- `N` is the FFT size (the number of points)

Each output bin combines contributions from every input sample, weighted by a complex exponential (a twiddle factor).

---

## Why This Kernel Matters

The FFT is one of the most widely used algorithms in signal processing, and it appears throughout real systems:

- Audio and speech processing, where signals are analyzed and filtered in the frequency domain
- Wireless communications (such as OFDM), where the FFT separates a signal into subcarriers
- Radar and sonar, where returns are analyzed for range and velocity
- Spectrum analysis and condition monitoring, where frequency content reveals system behavior
- Scientific computing and image processing, where transforms accelerate convolution and correlation

Because it combines complex arithmetic, non-sequential memory access, and multiple dependent stages, the FFT is a strong measure of how well an architecture handles structured, data-intensive workloads.

---

## Why Efficient Hardware Performs Well

The Electron E1x runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

This fits the radix-4 FFT well:

- Each radix-4 butterfly is a fixed pattern of complex multiplies and adds that maps cleanly onto the dataflow graph
- Intermediate results stay resident on the Fabric and pass directly from one butterfly to the next, so the transform runs largely in place
- The twiddle-factor and bit-reversal access patterns are known ahead of time, so addressing is baked into the graph rather than recomputed at runtime
- Complex multiply-accumulate operations run in parallel across the Fabric, keeping arithmetic throughput high
- The regular, predictable structure means very little energy goes to control overhead

The result is high signal-processing throughput at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## Configurable Parameters

These definitions in `main.c` control the benchmark. Change them to resize the problem or re-run it.

| Definition       | Default | Effect                                                                                                                                                                                                                                                          |
| ---------------- | ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `NUM_ITERATIONS` | `1`     | This is how many times the kernel runs. Increase it to average out measurement noise when benchmarking.                                                                                                                                                                 |
| `FFT_SIZE`       | `4096`  | This is the number of FFT points (`N`). This sets the problem size: a larger value means more input samples and more butterfly stages. If you change it, the build must be updated so the twiddle factors and expected reference data are regenerated for the new size. |

The input samples (`sample_input`) and the expected reference output (`expectedR`, `expectedI`) are generated for the chosen `FFT_SIZE`. Each output bin is compared against the reference with a tolerance of 10 per component to account for fixed-point quantization. If you change the size or the input data, the reference data must be regenerated to match the new correct result.
