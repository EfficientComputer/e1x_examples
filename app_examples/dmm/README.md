# Optimized Biquad IIR Filter

This example provides an **optimized Dense Matrix Multiply (DMM)** tailored for Efficient Computer (EFF) hardware.  
It is designed to showcase **real-world performance** and highlight the advantages of EFF’s architecture.

---

## 1. Overview

### What is a Dense Matrix Multiply?
A **dense matrix multiply** is a fundamental linear algebra function for the multiplication of dense matricies.

---

## 2. Why EFF Hardware Performs Well

This implementation is optimized to:

- Employ the maximum number of arithmetic units
- Parallelize computation and memory load/store 

---

## 3. Code Structure

### Core Files

- **`dmm.c`**
  - Reference DMM implementation and an EFF-optimized variant
  - `reformat(...)`: Reorders matrix A for efficient computation and memory access
  - `transpose_reformat(...)`: Transposes and reorders matrix B for efficient computation and memory access
  - `dmm_helper(...)`: Performs matrix multiplication using prepared matricies.
  - `dmm(...)`: DMM API used by the app
  - `dmm_reference(...)`: Reference DMM implementation useful for correctness checking

- **`main.c`**
  - Application entry point
  - Selects input size (via `EFFX_INPUT_SIZE`)
  - Generates sample matricies and runs the matrix multiplication
  - Intended for benchmarking + correctness checking

- **`CMakeLists.txt`**
  - Build configuration
  - Defines the `dmm` app and its build targets (sim / fabric / scalar_fabric)

---

## 4. Build Instructions

### Prerequisites
- CMake
- EFF SDK environment / toolchain

### Build Steps (typical EFF SDK flow)
Refer to the main build instruction

## 5. How to Use and Test

### Input data

Input size selection in `main.c` is controlled by `EFFX_INPUT_SIZE` (and defaults vary based on build mode).  
There are three size option:  64x64(small), 128x128(medium), 256x256(large)

### Running
Run the built binary using your normal EFF workflow (sim or fabric). The app:
1. Loads the input samples
2. Runs the matrix multiply using optimized and reference implementations
3. Checks optimized result against reference

### Performance benchmarking
This example is meant to highlight:
- instruction/memory efficiency of the optimized kernel

