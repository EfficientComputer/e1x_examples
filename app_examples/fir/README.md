# Finite Impulse Response Filter (FIR)

This example computes a **finite impulse response (FIR) filter** on the Electron E1 general-purpose processor. It is a compact, compute-dense signal-processing kernel used to show how a sliding-window dot product maps onto the Fabric architecture.

---

## 1. Overview

### What is a FIR filter?
A FIR filter produces each output sample as a weighted sum of the most recent input samples. It slides a fixed set of `W` filter coefficients (called taps) across the input signal, and at every position it multiplies each tap by the input sample under it and adds the products together. "Finite impulse response" means the output depends only on a finite window of past inputs, so the filter has no feedback and is always stable.

### Mathematical Definition

    o[i] = Σ w[j] * x[i + j]     for j = 0 .. W-1

Where:
- `x` is the input signal
- `w` is the vector of `W` filter taps (coefficients)
- `o` is the filtered output signal (`N` samples)

Each output sample `o[i]` costs `W` multiply-accumulate (MAC) operations, for `N * W` MACs total. Successive windows overlap heavily, so most of the input samples used for one output are reused for the next.

---

## 2. Why This Kernel Matters

FIR filtering is one of the most common operations in signal processing:

- **Audio processing**, where FIR filters implement equalizers, crossovers, and noise shaping
- **Communications and radio**, where they perform pulse shaping, channel filtering, and matched filtering
- **Sensor and biomedical front ends**, where they smooth, band-limit, and remove drift from sampled measurements
- **Image and video processing**, where separable FIR kernels blur, sharpen, and resample

Because it is arithmetic-bound and streams predictably over its input, FIR filtering is a good measure of how well an architecture sustains a tight multiply-accumulate loop with heavy data reuse.

---

## 3. Why EFF Hardware Performs Well

The Electron E1 runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

This fits FIR filtering well:

- The tap-by-tap multiply-accumulate operations for one output sample run **in parallel** across the Fabric
- The filter taps stay **resident on the Fabric** and are reused for every output sample instead of being reloaded
- The sliding window overlaps from one sample to the next, so each input value is **reused across many outputs** and memory traffic stays low
- The regular, predictable streaming access pattern maps cleanly onto a static dataflow graph, so very little energy goes to control overhead

The result is high arithmetic throughput at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## 4. Configurable Parameters

These definitions in `fir.h` and `main.c` control the benchmark. Change them to resize the problem or re-run it.

| Definition | Default | Effect |
|---|---|---|
| `N` | `64` | The number of output samples produced, defined in `fir.h`. A larger value means a longer output signal and more MACs per run. |
| `W` | `16` | The number of filter taps, defined in `fir.h`. It sets how many MACs go into each output sample and the length of the sliding window. |
| `NUM_ITERATIONS` | `1` | How many times the kernel runs. Increase it to average out measurement noise when benchmarking. |

The input signal `x` and the taps `w` are generated inline in `main.c` (`x[i] = i` and `w[i] = i - 8`), and the output is checked against a hardcoded expected value of `280` for `o[0]`. If you change `N`, `W`, or the way the inputs are generated, that expected value must be updated to match the new correct result.
