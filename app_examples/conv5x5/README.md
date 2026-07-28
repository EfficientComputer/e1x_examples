# 5x5 2D Convolution (Conv5x5)

This example computes a **5x5 two-dimensional convolution** of an input image with a fixed filter on the Electron E1x general-purpose processor. It is a compact, compute-dense image and signal-processing kernel used to show how a sliding-window stencil maps onto the Fabric architecture.

---

## Overview

### What is a 5x5 2D Convolution?

A 5x5 2D convolution slides a 5x5 filter (also called a kernel) across a 2D input image. At every output position it multiplies the 25 filter weights by the 25 input pixels under the window and sums the products into a single output pixel. The window then steps to the next position and the process repeats. This is the core operation behind image blurring, sharpening, edge detection, and the convolutional layers of neural networks.

### Mathematical Definition

    O(x,y) = Σ_i Σ_j I(x+i, y+j) * K(i,j)     for i,j = 0 .. 4

Where:

- `I` is the input image
- `K` is the 5x5 filter (25 weights)
- `O` is the output image
- `i` and `j` range over the 5x5 window

Each output pixel costs 25 multiply-accumulate (MAC) operations. For an output of `W` by `H` pixels this is `25 * W * H` MACs total.

---

## Why This Kernel Matters

2D convolution is one of the most common operations in real systems:

- **Image processing**, where filters blur, sharpen, denoise, and detect edges
- **Computer vision**, where convolutional neural networks extract features from images
- **Audio and signal processing**, where 2D filters operate on spectrograms and other time-frequency data
- **Scientific computing**, where stencil operations solve differential equations on grids

Because it applies the same small set of weights to a large stream of data, it is a good measure of how well an architecture reuses data and sustains arithmetic throughput while streaming an image through memory.

---

## Why Efficient Hardware Performs Well

The Electron E1x runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

This fits 5x5 2D convolution well:

- The 25 multiply-accumulate operations per output pixel run **in parallel** across the Fabric
- Filter weights stay **resident on the Fabric**, so they are loaded once and reused for every output pixel
- The window slides by reusing input pixels already held on the Fabric, so each pixel is **read from memory once** and shared across the outputs that overlap it
- Partial sums accumulate **on the Fabric** instead of spilling to memory
- Computation **overlaps with memory loads**, hiding load latency
- The regular, predictable sliding-window access pattern maps cleanly onto a static dataflow graph, so very little energy goes to control overhead

The result is high arithmetic throughput at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## Configurable Parameters

These definitions in `main.c` control the benchmark. Change them to resize the problem or re-run it.

| Definition         | Default | Effect                                                                                                                                                                |
| ------------------ | ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `NUM_ITERATIONS`   | `1`     | This is how many times the kernel runs. Increase it to average out measurement noise when benchmarking.                                                                       |
| `N`                | `64`    | This is the side length of the input image (both width and height). This sets the problem size: a larger value means a larger image and more MACs per run.                    |
| `N_PAD`            | `N + 4` | This is the padded side length of the input buffer. The filter is 5 wide, so the buffer is padded by 4 to hold the full sliding window. Derived from `N`; keep it as `N + 4`. |
| `RANDOMIZE_FILTER` | `1`     | When nonzero, the filter weights are overwritten with pseudo-random values at startup instead of the built-in weights.                                                |
| `RANGE`            | `10`    | This is the range of the pseudo-random filter weights. It also centers them around zero by subtracting `(RANGE - 1) / 2`.                                                     |
