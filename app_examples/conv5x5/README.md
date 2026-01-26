# Optimized 5×5 Convolution Kernel

This example provides an **optimized 5×5 convolution (Conv2D) implementation** tailored for **Efficient Computer (EFF) hardware**.  
It is designed to demonstrate **real-world ML and vision workloads** and to highlight the advantages of EFF’s architecture when executing **compute-dense and memory-intensive kernels**.

---

## 1. Overview

### What is a 5×5 Convolution?
A **5×5 convolution** is a fundamental operation in image processing and convolutional neural networks (CNNs).  
It computes each output pixel as a weighted sum of a **5×5 neighborhood** of input pixels.

Common use cases include:
- CNN feature extraction (early vision layers)
- Image filtering (blur, sharpen, edge detection)
- Embedded computer vision
- Preprocessing for ML inference

---

### Mathematical Definition

O(x, y) = Σ Σ I(x+i, y+j) * K(i, j)
i=0..4, j=0..4

Where:
- `I` is the input feature map
- `K` is the 5×5 convolution kernel
- `O` is the output feature map

Each output pixel requires **25 multiply-accumulate (MAC)** operations.

---

## 2. Why This Kernel Matters

5×5 convolution is a **realistic and challenging benchmark** because it:

- Has **high arithmetic intensity** (many MACs per output)
- Exhibits **sliding-window memory access**
- Is sensitive to **data layout and memory efficiency**
- Cannot be trivially vectorized without careful scheduling

These characteristics make it an excellent test of **architectural efficiency**, not just peak FLOPs.

---

## 3. Why EFF Hardware Performs Well

This implementation is optimized to:

- Keep kernel coefficients and partial sums in **registers**
- Use **tight inner loops** with predictable control flow
- Minimize redundant memory loads
- Exploit EFF’s ability to **overlap computation and memory operations**
- Sustain high throughput in nested loops

EFF hardware enables:
- Efficient handling of **loop-carried dependencies**
- Strong performance for **sliding-window access patterns**
- High utilization without heavy reliance on SIMD intrinsics

---

## 4. Code Structure

### Core Files

- **`conv5x5.c`**
  - Contains the reference and optimized 5×5 convolution implementations
  - Implements the core convolution kernel
  - Designed for clarity and performance comparison

- **`main.c`**
  - Application entry point
  - Initializes deterministic input data
  - Invokes both reference and optimized kernels
  - Compares outputs for correctness
  - Reports mismatches and performance information

---

## 5. Data Layout and Assumptions

- Input data is stored in **row-major order**
- Output data follows the same layout
- Kernel coefficients are fixed at compile time
- No padding

This layout mirrors typical CNN and image-processing pipelines.

---

## 6. Build Instructions

### Prerequisites
- CMake
- EFF SDK / toolchain
- Standard EFF build environment (sim / fabric)

### Build
Refer to the main EFF SDK build instructions.  
The build process compiles:
- Reference convolution
- Optimized convolution
- Test harness in `main.c`

---

## 7. How to Run and Test

### Running
Run the built binary using your normal EFF workflow (sim or fabric).

The application performs the following steps:
1. Initializes input image data
2. Runs the reference 5×5 convolution
3. Runs the optimized 5×5 convolution
4. Compares output buffers
5. Reports mismatches (if any)

---

### Correctness Checking

Correctness is validated by comparing each output element:
ref_out[] vs test_out[]

---

## 8. Performance Benchmarking

This example is intended to highlight:
- Sustained throughput for **sliding-window convolution**
- Instruction-level efficiency in nested loops
- Memory efficiency under high data reuse
- EFF’s advantages for **ML and vision kernels**

Performance can be measured using:
- Cycle counters
- EFF simulator statistics
- Hardware profiling tools

---

## 9. Why This Example Is Useful

5×5 convolution is representative of:
- Early layers in CNNs (e.g., MobileNet-style models)
- Embedded vision workloads
- Image preprocessing pipelines

Efficient execution of this kernel directly impacts:
- End-to-end ML inference latency
- Power consumption
- Real-time system responsiveness

---

## 10. Summary

- Demonstrates an optimized **5×5 convolution kernel**
- Validates correctness against a reference implementation
- Highlights **compute and memory efficiency** on EFF hardware
- Relevant to ML, vision, and DSP workloads
- Shows why EFF excels at real-world kernels, not just microbenchmarks