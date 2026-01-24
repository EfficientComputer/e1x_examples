# JPEG Decode and Reference Validation Example

This code package provides a **JPEG decoding example** designed for **Efficient Computer (EFF) hardware and toolchains**.  
It demonstrates **end-to-end JPEG decoding**, reference validation, and correctness checking using a **known-good software reference decoder**.

The example is intended to showcase:
- Realistic image-processing workloads
- Memory and compute behavior of JPEG decode
- Correctness validation against a trusted reference

---

## 1. Overview

### What Does This Example Do?
This application:
1. Decodes a JPEG image using an optimized JPEG decoder
2. Decodes the same JPEG using a **reference implementation**
3. Compares the decoded RGB output buffers
4. Reports mismatches and correctness status

This mirrors a **real-world validation workflow** used when porting or optimizing codecs for new hardware.

---

## 2. Why JPEG Decode Is Important

JPEG decoding is a **widely deployed, real-world workload** found in:
- Cameras and imaging pipelines
- Embedded vision systems
- ML preprocessing pipelines
- Mobile and edge devices

It is challenging because it combines:
- Bitstream parsing
- Huffman decoding
- IDCT and dequantization
- Color conversion
- Non-trivial memory access patterns

As such, JPEG decode is an excellent benchmark for **compute efficiency, memory bandwidth, and control flow handling**.

---

## 3. Code Structure

### Core Files

- **`jpeg.c`**
  - Optimized JPEG decoder implementation
  - Handles JPEG bitstream parsing and image reconstruction
  - Produces decoded image output in RGB format

- **`stbi_ref.h`**
  - Reference JPEG decoder based on `stb_image`
  - Used as a **golden reference** for correctness
  - Not optimized; intended only for validation

- **`image.c`**
  - Contains embedded JPEG image data
  - Defines image size, format, and input buffer
  - Used as deterministic test input

- **`generate_answer.c`**
  - Runs the reference decoder
  - Generates expected (golden) output
  - Can be used offline or as part of the validation flow

- **`main.c`**
  - Application entry point
  - Invokes both optimized and reference decoders
  - Compares decoded RGB outputs
  - Reports mismatches and summary results

---

## 4. Data Flow

JPEG bitstream (image.h)
|
v
Optimized JPEG Decoder (jpeg.c)
|
v
Decoded RGB output (opt)

JPEG bitstream (image.h)
|
v
Reference Decoder (stbi_ref.h)
|
v
Decoded RGB output (ref)

opt vs ref  → correctness check

---

## 5. Output Format

- Decoded output is **RGB888**
- Data layout: R0, G0, B0, R1, G1, B1, …

- Output buffers are compared **element-by-element**
- Floating-point or integer comparisons use appropriate tolerances (if applicable)

---

## 6. Correctness Validation

Correctness is verified by comparing the optimized decoder output against the reference decoder:

- Byte-wise comparison for RGB data
- Detailed mismatch reporting:
- Pixel index
- Channel
- Expected vs actual value

This ensures that:
- Bitstream parsing is correct
- IDCT and color conversion are accurate
- No memory corruption or alignment issues exist

---

## 7. Build Instructions

### Prerequisites
- C compiler compatible with EFF SDK
- Standard EFF build environment (sim / fabric)
- CMake or existing EFF build system

### Build
Refer to the standard EFF SDK build instructions.  
The build compiles:
- Optimized JPEG decoder
- Reference decoder
- Test harness (`main.c`)

---

## 8. Running the Example

Run the compiled binary using your standard EFF workflow.

At runtime, the program:
1. Loads the embedded JPEG image
2. Decodes using the optimized JPEG path
3. Decodes using the reference path
4. Compares outputs
5. Prints pass/fail status and mismatch details (if any)

---

## 9. Performance Considerations

This example highlights:
- Control-heavy workloads with irregular memory access
- Sustained throughput of IDCT and color conversion
- Flash-to-SRAM data movement costs
- Benefits of optimized decode kernels on EFF hardware

Performance can be evaluated using:
- Cycle counters
- Simulator statistics
- Hardware profiling tools

---

## 10. Why This Example Matters

JPEG decoding represents:
- A **real customer workload**
- A mix of scalar, control, and arithmetic operations
- A stress test for memory systems and instruction scheduling

Efficient execution directly impacts:
- Image decode latency
- Power consumption
- End-to-end vision pipeline performance

---

## 11. Summary

- Demonstrates a full **JPEG decode pipeline**
- Validates correctness against a trusted reference
- Uses deterministic embedded test data
- Highlights EFF hardware strengths on real-world codecs
- Suitable for benchmarking, validation, and regression testing

---

## 12. Notes

- Reference decoder is included for validation only
- Optimized decoder behavior should always match reference output
- Any mismatch indicates a functional bug and should be investigated


## 13. generate_answer.c

- To build:  "gcc -O3 -Wall generate_answer.c image.c -o generate_answer". 
- To run:    "./generate_answer" 
- Output :   answer.jpg — generated JPEG image
             answer.c — auto-generated C source with embedded JPEG data

---