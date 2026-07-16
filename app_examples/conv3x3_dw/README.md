# Depthwise 3x3 Convolution (Conv3x3-DW)

This example computes a **depthwise 3x3 convolution** on the Electron E1x general-purpose processor. It applies one independent 3x3 filter to each input channel, a pattern that dominates the compute of MobileNet-style convolutional networks.

---

## 1. Overview

### What is a depthwise 3x3 convolution?

A standard convolution layer mixes information across channels. A depthwise convolution does not: it keeps each input channel separate and applies its own 3x3 filter to that channel alone. With `C` channels there are `C` independent 3x3 filters, and channel `c` of the output depends only on channel `c` of the input. This makes depthwise convolution far cheaper than a full multi-channel convolution, which is why it is a standard building block in efficient mobile vision models.

This example uses int8 inputs and int8 filter weights, accumulates into int32, uses a stride of 1, and applies SAME padding (a one-element border of zeros on every side) so the output has the same height and width as the input.

### Mathematical Definition

    O(x, y, c) = Σ_i Σ_j I(x + i, y + j, c) * K(i, j, c)     for i, j = 0 .. 2

Where:

- `I` is the input, with `C` channels
- `K(i, j, c)` is the 3x3 filter for channel `c`
- `O` is the output, with `C` channels
- each output channel `c` uses only input channel `c`

Each output value costs nine multiply-accumulate (MAC) operations, and there are `C` independent channels.

---

## 2. Why This Kernel Matters

Depthwise 3x3 convolution is a workhorse of modern efficient neural networks:

- **Mobile and edge vision models**, where MobileNet and similar architectures use depthwise convolution to cut cost dramatically
- **On-device inference**, where per-channel spatial filtering runs within a tight energy budget
- **Separable convolution designs**, where a depthwise stage does spatial mixing and a following pointwise stage mixes channels
- **Real-time image pipelines**, where lightweight per-channel filters process camera frames continuously

Because each channel is independent and the same small filter is reused across every spatial position, it is a good measure of how well an architecture exploits parallelism and data reuse under low precision.

---

## 3. Why EFF Hardware Performs Well

The Electron E1x runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

This fits depthwise 3x3 convolution well:

- The channels are **independent**, so their work can proceed in parallel across the Fabric
- The nine multiply-accumulate operations per output run in parallel, and the int8 inputs keep data movement small
- Per-channel filter taps stay **resident on the Fabric** and are reused for every output position in that channel
- As the window slides sideways, previously loaded values are **reused** instead of reloaded, so only the new column is fetched
- Computation **overlaps with memory loads**, hiding load latency
- The regular, predictable access pattern maps cleanly onto a static dataflow graph, so very little energy goes to control overhead

The result is high arithmetic throughput at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## 4. Configurable Parameters

These definitions in `main.c` and the app's `conv3x3_dw.h` header control the benchmark. Change them to resize the problem or re-run it. The optimized kernel is compared against a reference implementation, so correctness holds automatically when you resize.

| Definition       | Default | Effect                                                                                          |
| ---------------- | ------- | ----------------------------------------------------------------------------------------------- |
| `NUM_ITERATIONS` | `1`     | How many times the kernel runs. Increase it to average out measurement noise when benchmarking. |
| `M`              | `48`    | The number of input rows (image height).                                                        |
| `N`              | `48`    | The number of input columns (image width).                                                      |
| `NCHANNELS`      | `8`     | The number of channels. Each channel gets its own independent 3x3 filter.                       |
| `KERNEL_DIM`     | `3`     | The filter side length. This kernel is specialized for 3x3 filters.                             |
| `STRIDE_H`       | `1`     | The vertical stride between output positions.                                                   |
| `STRIDE_W`       | `1`     | The horizontal stride between output positions.                                                 |
| `PAD_TOP`        | `1`     | Rows of zero padding added at the top.                                                          |
| `PAD_LEFT`       | `1`     | Columns of zero padding added on the left.                                                      |
| `PAD_RIGHT`      | `1`     | Columns of zero padding added on the right.                                                     |
| `PAD_BOTTOM`     | `1`     | Rows of zero padding added at the bottom.                                                       |
| `OH`             | `M`     | The output height. With SAME padding and stride 1 it equals the input height.                   |
| `OW`             | `N`     | The output width. With SAME padding and stride 1 it equals the input width.                     |
| `OUT_MUL`        | `1`     | The output multiplier per channel (one output channel per input channel for depthwise).         |
