# Multi-Channel 3x3 Convolution (Conv3x3xN)

This example computes a **multi-channel 3x3 convolution** on the Electron E1 general-purpose processor. It convolves a 3x3 filter over every input channel and sums the results across channels, which is a full convolutional layer.

---

## 1. Overview

### What is a multi-channel 3x3 convolution?
A multi-channel (N-channel) convolution takes an input with several channels and, at each output position, convolves a 3x3 filter over all of them and adds the results together into a single output value. Unlike a depthwise convolution, which keeps channels separate, this operation mixes information across channels: every output value combines contributions from all `N` input channels. This is the standard convolutional layer used throughout vision models.

This example uses int32 data and slides the 3x3 window across the image, accumulating a running sum over the channels for each output position.

### Mathematical Definition

    O(x, y) = Σ_k Σ_i Σ_j I(x + i, y + j, k) * K(i, j, k)     for i, j = 0 .. 2, k = 0 .. N-1

Where:
- `I` is the input, with `N` channels
- `K(i, j, k)` is the 3x3 filter for channel `k`
- `O` is the output, summed across all channels

Each output value costs nine multiply-accumulate (MAC) operations per channel, for `9 * N` MACs total.

---

## 2. Why This Kernel Matters

Multi-channel 3x3 convolution is the core of most convolutional neural networks:

- **Convolutional neural networks**, where stacked multi-channel 3x3 layers do the bulk of feature extraction in vision models
- **Image classification and detection**, where each layer combines many input feature maps into new ones
- **Video and sensor fusion**, where multiple input planes are convolved and merged
- **Feature extraction pipelines**, where cross-channel mixing builds richer representations layer by layer

Because it combines dense arithmetic, cross-channel accumulation, and heavy window reuse, it is a good measure of how well an architecture sustains throughput while moving data efficiently.

---

## 3. Why EFF Hardware Performs Well

The Electron E1 runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

This fits multi-channel 3x3 convolution well:

- The many multiply-accumulate operations across the 3x3 window and the channels run **in parallel** across the Fabric
- Filter weights stay **resident on the Fabric** and are reused for every output position
- As the window slides sideways, previously loaded values are **reused** instead of reloaded, so only the new column is fetched
- Cross-channel partial sums accumulate on the Fabric while loads for the next window **overlap with computation**, hiding load latency
- The regular, predictable access pattern maps cleanly onto a static dataflow graph, so very little energy goes to control overhead

The result is high arithmetic throughput at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## 4. Configurable Parameters

These definitions in `main.c` and the app's `conv3x3xn.h` header control the benchmark. Change them to resize the problem or re-run it. The optimized kernel is compared against a reference implementation, so correctness holds automatically when you resize.

| Definition | Default | Effect |
|---|---|---|
| `NUM_ITERATIONS` | `1` | How many times the kernel runs. Increase it to average out measurement noise when benchmarking. |
| `M` | `16` | The number of input rows (image height). Keeping it a power of two helps performance. |
| `N` | `16` | The number of input columns (image width). Keeping it a power of two helps performance. |
| `KERNEL_DIM` | `3` | The filter side length. This kernel is specialized for 3x3 filters. |
| `NCHANNELS` | `4` | The number of input channels summed together for each output value. |
| `FASTPATH` | `1` | When enabled, adds padding columns to the input and output buffers so the store does not need to be gated, which improves performance. |
| `INSTRIDE` | `N + 1` | The row stride of the input buffer in elements. With `FASTPATH` enabled it includes one padding column to avoid memory bank conflicts. |
| `OUTSTRIDE` | `N + KERNEL_DIM - 1` | The row stride of the output buffer in elements. With `FASTPATH` enabled it includes the leading padding columns. |
