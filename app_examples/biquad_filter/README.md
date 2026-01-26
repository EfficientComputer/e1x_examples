# Optimized Biquad IIR Filter

This example provides an **optimized Biquad IIR filter implementation** tailored for Efficient Computer (EFF) hardware.  
It is designed to showcase **real-world DSP performance** and highlight the advantages of EFF’s architecture compared to traditional CPU, DSP, and SIMD-heavy solutions.

---

## 1. Overview

### What is a Biquad IIR Filter?
A **biquad filter** is a second-order Infinite Impulse Response (IIR) digital filter widely used in:

- Audio processing (EQ, filters)
- Sensor signal conditioning
- Communications
- Control systems

It is called *biquad* because it implements **two poles and two zeros**.

### Difference Equation
The filter follows the standard biquad equation:

y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2]
     - a1*y[n-1] - a2*y[n-2]


Where:
- `x[n]` is the input signal
- `y[n]` is the output signal
- `b0, b1, b2` are feed-forward coefficients (zeros)
- `a1, a2` are feedback coefficients (poles)

### Why This Matters
Unlike FIR filters, IIR filters have **loop-carried dependencies**, meaning:
- Each output depends on previous outputs
- Aggressive vectorization is limited or impossible

This makes IIR filters an excellent **real-world benchmark** for evaluating architectural efficiency.

---

## 2. Why EFF Hardware Performs Well

This implementation is optimized to:

- Keep filter state in **registers**
- Use tight loops with predictable control flow
- Exploit EFF’s ability to **pipeline computation** efficiently even with loop-carried dependencies
- Parallelize computation and memory load/store 


---

## 3. Code Structure

### Core Files

- **`biquad_filter.c`**
  - Reference biquad implementation and an EFF-optimized variant
  - `biquad_filter_q15(...)`: Q15 fixed-point kernel (reference vs optimized via compile-time switch)
  - `biquad_filter(...)`: wrapper API used by the app

- **`main.c`**
  - Application entry point
  - Selects input size (e.g., via `EFFX_INPUT_SIZE`)
  - Loads deterministic sample data and runs the filter
  - Intended for benchmarking + correctness checking

- **`gen_data.py`**
  - Generates deterministic input samples as a C include file
  - Produces `data_<N>.h.inc` given `--op-count <N>`

- **`CMakeLists.txt`**
  - Build configuration + data-generation hooks
  - Defines the `biquad_4samples` app and its build targets (sim / fabric / scalar_fabric)
  - Integrates `gen_data.py` to generate input headers for multiple benchmark sizes

---

## 4. Build Instructions

### Prerequisites
- Python 3
- CMake
- EFF SDK environment / toolchain (so that `add_eff_app(...)` targets are available)

### What the build does
The build first generates sample input headers via `gen_data.py` and then compiles the app.  
`CMakeLists.txt` generates multiple datasets (small/medium/large) and wires them into build targets.

### Build Steps (typical EFF SDK flow)
Refer to the main build instruction


## 5. How to Use and Test

### Input data
The program includes a generated header for input samples:
- `biquad_filter_input.h.inc`

Input size selection in `main.c` is controlled by `EFFX_INPUT_SIZE` (and defaults vary based on build mode).  
There are three size option:  512(small), 4096(medium), 16384(large)

### Running
Run the built binary using your normal EFF workflow (sim or fabric). The app:
1. Loads the input samples
2. Runs the biquad filter
3. Writes results to the `output[]` buffer


### Correctness
To validate correctness:
- Run the same input through the reference implementation and compare output samples.
- For fixed-point variants, allow for expected quantization differences depending on coefficient scaling / rounding.

### Performance benchmarking
This example is meant to highlight:
- sustained throughput with loop-carried dependencies
- instruction/memory efficiency of the optimized kernel
- EFF HW achieves 7.3 cycles per sample per stage



