# 3x3 Convolution (Conv3x3)

This example computes a **single-channel 3x3 2D convolution** on the Electron E1 general-purpose processor. It slides a 3x3 filter over a 2D input image and produces one output value per position, a core building block of image and signal processing.

---

## 1. Overview

### What is a 3x3 convolution?
A 3x3 convolution slides a small 3x3 filter (also called a kernel) across a 2D input image. At every position, it multiplies the nine overlapping input values by the nine filter weights and sums the results into a single output value. The filter then shifts by one column and repeats, so neighboring outputs share most of their input values. This example is single-channel: it works on one 2D plane of data with one 3x3 filter.

### Mathematical Definition

    O(x, y) = Σ_i Σ_j I(x + i, y + j) * K(i, j)     for i, j = 0 .. 2

Where:
- `I` is the input image
- `K` is the 3x3 filter (nine weights)
- `O` is the output image

Each output value costs nine multiply-accumulate (MAC) operations, one per filter weight.

---

## 2. Why This Kernel Matters

3x3 convolution is one of the most widely used operations in image and signal processing:

- **Image filtering**, where a 3x3 kernel blurs, sharpens, or denoises a picture
- **Edge and feature detection**, where Sobel and similar 3x3 operators highlight gradients and boundaries
- **Convolutional neural networks**, where stacked 3x3 layers are the dominant compute cost in vision models
- **Computer vision on sensors**, where local 3x3 filters run close to the camera before higher-level processing

Because the same small filter is applied everywhere and neighboring windows overlap heavily, it is a good measure of how well an architecture reuses data and sustains dense arithmetic.

---

## 3. Why EFF Hardware Performs Well

The Electron E1 runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

This fits 3x3 convolution well:

- The nine multiply-accumulate operations per output run **in parallel** across the Fabric
- Filter weights stay **resident on the Fabric** and are reused for every output position
- As the window slides sideways, previously loaded values are **reused** instead of reloaded, so only the new column is fetched
- Computation **overlaps with memory loads**, hiding load latency
- The regular, predictable access pattern maps cleanly onto a static dataflow graph, so very little energy goes to control overhead

The result is high arithmetic throughput at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## 4. Configurable Parameters

These definitions in `main.c` and the app's `conv3x3.h` header control the benchmark. Change them to resize the problem or re-run it.

| Definition | Default | Effect |
|---|---|---|
| `NUM_ITERATIONS` | `1` | How many times the kernel runs. Increase it to average out measurement noise when benchmarking. |
| `M` | `64` | The number of input rows (image height). Larger values mean a larger image and more output positions. |
| `N` | `64` | The number of input columns (image width). Larger values mean a larger image and more output positions. |
| `INSTRIDE` | `N` | The row stride of the input buffer in elements. It sets how far apart successive input rows are in memory. |
| `OUTSTRIDE` | `N` | The row stride of the output buffer in elements. It sets how far apart successive output rows are in memory. |
