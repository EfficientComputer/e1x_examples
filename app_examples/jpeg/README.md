# Optimized JPEG Encoder

This example provides an **optimized JPEG image encoder** tailored for Efficient Computer (EFF) hardware.  
It is designed to showcase **real-world image compression performance** and highlight the advantages of EFF's architecture for compute-intensive multimedia workloads.

---

## 1. Overview

### What is JPEG Encoding?
**JPEG (Joint Photographic Experts Group)** is one of the most widely used image compression formats. The encoding process transforms raw pixel data into a compressed bitstream using a combination of:

- Color space conversion (RGB to YCbCr)
- Chroma subsampling
- Discrete Cosine Transform (DCT)
- Quantization
- Huffman encoding

Common use cases include:
- Digital cameras and smartphones
- Web image serving
- Medical imaging
- Embedded vision systems
- IoT edge devices

---

### Encoding Pipeline

1. **Color Conversion**: RGB → YCbCr (luminance + chrominance)
2. **Block Division**: Image split into 8×8 or 16×16 macroblocks
3. **DCT**: 2D Discrete Cosine Transform on each block
4. **Quantization**: Frequency coefficients divided by quality-dependent tables
5. **Entropy Coding**: Huffman encoding of quantized coefficients

---

## 2. Why This Kernel Matters

JPEG encoding is a **realistic and challenging benchmark** because it:

- Involves **multiple processing stages** with different computational characteristics
- Requires **color space transformations** with floating-point or fixed-point arithmetic
- Uses **2D DCT** which is computationally intensive
- Has **variable-length encoding** with bit manipulation
- Exhibits **block-based processing** with good data locality

These characteristics make it an excellent test of **end-to-end system efficiency**, not just peak compute.

---

## 3. Why EFF Hardware Performs Well

This implementation is optimized to:

- Use **fixed-point arithmetic** (Q16 format) for efficient integer computation
- Process **multiple macroblocks in batches** to amortize setup overhead
- Exploit EFF's ability to **parallelize nested loops** in DCT and color conversion
- Keep quantization tables and DCT coefficients in **registers**
- Use **tight inner loops** for Huffman encoding
- Minimize memory bandwidth with **in-place processing**

EFF hardware enables:
- Efficient execution of **multi-stage pipelines**
- High throughput for **2D transforms**
- Strong performance for **bit manipulation** in entropy coding

---

## 4. Code Structure

### Core Files

- **`jpeg.c`**
  - Complete JPEG encoder implementation (based on stb_image_write)
  - `stbiw__jpg_DCT(...)`: 1D DCT butterfly computation
  - `stbiw__convert_colors_16/8(...)`: RGB to YCbCr color conversion
  - `stbiw__subsample_UV(...)`: Chroma subsampling (4:2:0)
  - `stbiw__jpg_process_dct_rows/cols(...)`: 2D DCT via row-column decomposition
  - `stbiw__jpg_encode_DC_AC(...)`: Huffman encoding of DC and AC coefficients
  - `stbi_write_jpg(...)`: Main encoder API

- **`image.c` / `image.h`**
  - Test image data (embedded as C array)
  - Image dimensions and format definitions

- **`main.c`**
  - Application entry point
  - Runs the JPEG encoder (10 iterations for profiling)
  - Validates output size against expected reference
  - Optionally writes output to file in simulation mode

- **`CMakeLists.txt`**
  - Build configuration
  - Defines the `jpeg` app and its build targets (sim / fabric / scalar)
  - Supports both E1 and E1x architectures

---

## 5. Specification

- **Input Format**: Raw RGB image (3 components per pixel)
- **Output Format**: Standard JPEG bitstream
- **Quality Setting**: 50 (configurable, enables 4:2:0 chroma subsampling)
- **Arithmetic**: Q16 fixed-point
- **Accuracy**: Output size within 2% of reference (accounts for float vs fixed-point differences)

---

## 6. Build Instructions

### Prerequisites
- CMake
- EFF SDK environment / toolchain

### Build Steps (typical EFF SDK flow)
Refer to the main build instructions.

The build compiles:
- Optimized JPEG encoder
- Test image data
- Correctness validation in `main.c`

---

## 7. How to Run and Test

### Input Data
The program uses an embedded test image:
- Defined in `image.c` / `image.h`
- RGB format, 3 bytes per pixel

### Running
Run the built binary using your normal EFF workflow (sim or fabric). The app:
1. Loads the source image
2. Runs the JPEG encoder (10 iterations for profiling)
3. Validates output size against reference
4. Reports PASS/FAIL status
5. (Simulation only) Writes `test.jpg` to disk for visual inspection

### Correctness Checking
Correctness is validated by comparing output file size:

```
abs(output_arr_length - answer_length) <= answer_length / 50
```

A 2% tolerance accounts for minor differences between floating-point and fixed-point implementations.

### Visual Verification
In simulation mode, the encoder writes `test.jpg` which can be opened with any image viewer to verify visual quality.

---

## 8. Performance Benchmarking

This example is intended to highlight:
- End-to-end **image compression throughput**
- Efficiency of **multi-stage processing pipelines**
- Performance of **2D DCT computation**
- Bit manipulation efficiency in **Huffman encoding**
- EFF's advantages for **multimedia workloads**

Performance can be measured using:
- Cycle counters
- EFF simulator statistics
- Hardware profiling tools

---

## 9. Why This Example Is Useful

JPEG encoding is representative of:
- Camera and imaging pipelines
- Video codec intra-frame compression
- Web and cloud image processing
- Embedded vision applications

Efficient execution of JPEG encoding directly impacts:
- Camera capture latency
- Power consumption in mobile devices
- Throughput in image processing servers
- Real-time video streaming quality

---

## 10. Summary

- Demonstrates an optimized **JPEG encoder**
- Uses **fixed-point arithmetic** for efficient computation
- Implements complete pipeline: color conversion, DCT, quantization, Huffman coding
- Validates correctness against reference output
- Highlights **multimedia processing efficiency** on EFF hardware
- Relevant to imaging, video, and embedded vision workloads
