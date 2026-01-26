# Optimized FFT (Fast Fourier Transform)

This example provides an **optimized radix-4 FFT implementation** tailored for Efficient Computer (EFF) hardware.  
It is designed to showcase **real-world DSP performance** and highlight the advantages of EFF's architecture for frequency-domain signal processing.

---

## 1. Overview

### What is an FFT?
The **Fast Fourier Transform (FFT)** is an efficient algorithm for computing the Discrete Fourier Transform (DFT), converting signals between the time domain and the frequency domain.

Common use cases include:
- Audio and speech processing
- Spectrum analysis
- Communications (OFDM, channel estimation)
- Radar and sonar signal processing
- Image processing
- Scientific computing

---

### Mathematical Definition

The DFT of a sequence x[n] of length N is:

X[k] = Σ x[n] * e^(-j*2π*k*n/N), for n = 0 to N-1

Where:
- `x[n]` is the input signal in time domain
- `X[k]` is the output in frequency domain
- `k` is the frequency bin index
- `N` is the FFT size (4096 in this example)

The FFT reduces computational complexity from O(N²) to O(N log N) by exploiting symmetry in the twiddle factors.

---

## 2. Why This Kernel Matters

FFT is a **fundamental DSP building block** and presents several optimization challenges:

- **Complex arithmetic**: Each butterfly operation involves complex multiplications
- **Bit-reversed addressing**: Output ordering requires non-sequential memory access
- **Twiddle factor lookups**: Pre-computed coefficients must be efficiently accessed
- **Multi-stage computation**: Radix-4 implementation requires log₄(N) stages

This makes FFT an excellent benchmark for evaluating **architectural efficiency** in DSP workloads.

---

## 3. Why EFF Hardware Performs Well

This implementation is optimized to:

- Use a **radix-4 butterfly** structure for fewer stages and better arithmetic density
- Keep twiddle factors and intermediate values in **registers**
- Exploit EFF's ability to **pipeline computation** across butterfly stages
- Use **Q15 fixed-point arithmetic** for efficient integer processing
- Implement **bit-reversal** using hardware-supported bit manipulation intrinsics
- Minimize memory bandwidth with **in-place computation**

EFF hardware enables:
- Efficient handling of **complex multiply-accumulate operations**
- High throughput for **multi-stage iterative algorithms**
- Strong performance without heavy reliance on SIMD intrinsics

---

## 4. Code Structure

### Core Files

- **`fft.c`**
  - Contains the optimized radix-4 FFT implementation
  - `fft_init_dst(...)`: Performs bit-reversal reordering of input data
  - `kiss_fft_run_layer(...)`: Executes a single stage of radix-4 butterflies
  - `fft4(...)`: Main FFT API that orchestrates all stages

- **`fft.h`**
  - FFT data type definitions (`fft_cpx` for complex numbers)
  - Function declarations
  - Sample input data (when `DEFINE_SAMPLE_INPUT` is set)

- **`fft_twiddles.h`**
  - Pre-computed twiddle factors for the 4096-point FFT
  - Twiddle schedule for each FFT stage

- **`expected.c`**
  - Reference output values for correctness validation
  - Generated from the unmodified KissFFT library

- **`main.c`**
  - Application entry point
  - Runs the FFT kernel (10 iterations for profiling)
  - Validates output against expected reference values

- **`CMakeLists.txt`**
  - Build configuration
  - Defines the `fft` app and its build targets (sim / fabric / scalar)

---

## 5. Specification

- **FFT Size**: 4096 bins (complex)
- **Data Format**: 16-bit Q15 fixed-point integers
- **Output Ordering**: Bit-reversed
- **Accuracy**: Element error must be < 10 compared to reference (to account for fixed-point quantization)

---

## 6. Build Instructions

### Prerequisites
- CMake
- EFF SDK environment / toolchain

### Build Steps (typical EFF SDK flow)
Refer to the main build instructions.

The build compiles:
- Optimized radix-4 FFT kernel
- Correctness validation in `main.c`

---

## 7. How to Run and Test

### Input Data
The program uses deterministic sample input defined in `fft.h`:
- Complex input array in Q15 format
- Each element represents a PCM sample

### Running
Run the built binary using your normal EFF workflow (sim or fabric). The app:
1. Loads the input samples
2. Runs the radix-4 FFT (10 iterations for profiling)
3. Compares output against reference values
4. Reports PASS/FAIL status

### Correctness Checking
Correctness is validated by comparing each output bin:

```
expectedR[i] vs out_buf[i].r  (real component)
expectedI[i] vs out_buf[i].i  (imaginary component)
```

The reference values are generated using the unmodified KissFFT library.

---

## 8. Performance Benchmarking

This example is intended to highlight:
- Sustained throughput for **multi-stage FFT computation**
- Efficient **complex arithmetic** handling
- Memory efficiency with **bit-reversed addressing**
- EFF's advantages for **DSP kernels**

Performance can be measured using:
- Cycle counters
- EFF simulator statistics
- Hardware profiling tools

---

## 9. Why This Example Is Useful

FFT is foundational for:
- Real-time audio processing
- Wireless communications (5G, WiFi, LTE)
- Embedded spectrum analyzers
- Edge AI preprocessing

Efficient execution of FFT directly impacts:
- End-to-end signal processing latency
- Power consumption in battery-powered devices
- Real-time system responsiveness

---

## 10. Summary

- Demonstrates an optimized **radix-4 FFT kernel** (4096-point)
- Uses **Q15 fixed-point arithmetic** for efficient computation
- Validates correctness against a reference implementation
- Highlights **DSP efficiency** on EFF hardware
- Relevant to audio, communications, radar, and scientific computing workloads
