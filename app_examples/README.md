# E1x App Examples

## Overview

This project provides example apps demonstrating optimized algorithms running on the E1x processor using the E1x Evaluation Kit (EVK). These examples serve as a reference for developers looking to understand how to leverage E1x's architecture for real-world DSP, image processing, communications, and compute workloads. Each app includes complete source code showing EFF-optimized implementations along with correctness validation.

---

## Apps

- **bfs** - Breadth-First Search (Graph)
- **biquad_filter** - Optimized Biquad IIR Filter (DSP)
- **cholesky_decomp** - Cholesky Decomposition (Linear Algebra)
- **conv3x3** - 3×3 Convolution Kernel (ML/Vision)
- **conv3x3_dw** - 3×3 Depthwise Convolution (ML/Vision)
- **conv3x3_dw_vec** - Vectorized 3×3 Depthwise Convolution (ML/Vision)
- **conv3x3_vec** - Vectorized 3×3 Convolution Kernel (ML/Vision)
- **conv3x3xn** - Multi-Channel 3×3 Convolution (ML/Vision)
- **conv3x3xn_vec** - Vectorized Multi-Channel 3×3 Convolution (ML/Vision)
- **conv5x5** - Optimized 5×5 Convolution Kernel (ML/Vision)
- **conv5x5_vec** - Vectorized 5×5 Convolution Kernel (ML/Vision)
- **conv5x5sym** - Symmetric 5×5 Convolution Kernel (ML/Vision)
- **conv5x5sym_vec** - Vectorized Symmetric 5×5 Convolution Kernel (ML/Vision)
- **crc32** - CRC-32 Checksum (Data Integrity)
- **dijkstras** - Dijkstra's Shortest Path (Graph)
- **dmm** - Dense Matrix Multiply (Linear Algebra)
- **dmm_i8** - Int8 Dense Matrix Multiply (Linear Algebra)
- **dmm_vectorized** - Vectorized Dense Matrix Multiply (Linear Algebra)
- **dmv** - Dense Matrix-Vector Multiply (Linear Algebra)
- **dmv_fp** - Floating-Point Dense Matrix-Vector Multiply (Linear Algebra)
- **dmv_vec** - Vectorized Dense Matrix-Vector Multiply (Linear Algebra)
- **fft** - Fast Fourier Transform (DSP)
- **fir** - FIR Filter (DSP)
- **fir_vec** - Vectorized FIR Filter (DSP)
- **hello** - Minimal "Hello World" example
- **isqrt** - Integer Square Root (Math)
- **jpeg** - JPEG Image Encoder (Multimedia)
- **knapsack** - 0/1 Knapsack Solver (Dynamic Programming)
- **ldpc** - LDPC Encoder/Decoder (Communications)
- **levenshtein** - Levenshtein Edit Distance (Dynamic Programming)
- **qr_decomp** - QR Decomposition (Linear Algebra)
- **quickstart** - Simple "Hello World" example
- **recursive_factorial** - Recursive Factorial (Math)
- **smv** - Sparse Matrix-Vector Multiply (Sparse Linear Algebra)
- **spadd** - Sparse Matrix Addition (Sparse Linear Algebra)
- **sparse** - Sparse Vector Operations (Sparse Linear Algebra)
- **spmspv** - Sparse Matrix–Sparse Vector Multiply (Sparse Linear Algebra)
- **sqrt_newton** - Newton's Method Square Root (Math)

## Getting Started

If you haven't already, please make sure you've set up your board and development environment using our [Evaluation Kit Setup Instructions](https://docs.efficient.computer/evaluation-kit). This project supports both our physical EVK and our Cloud EVK.

Additional details are contained within each app's folder.

### App Configuration

Each app is contained in its own directory with a `CMakeLists.txt` for building. Refer to the individual app README and source files for detailed usage and configuration options.

Many apps ship with both a portable reference implementation and a hand-optimized version tuned for the E1x Fabric. To build the hand-optimized versions, edit the top-level `CMakeLists.txt` in the `app_examples` folder and uncomment the following line:

```
add_compile_definitions(EFF_BLD_HAND_OPTIMIZED)
```

Then rebuild, as directed below.

### Building

To build all apps, execute the following commands from the top level of this folder:

```
mkdir bld
cd bld
cmake -G Ninja .. -DEFF_STDIO_PORT=3
ninja
```

This will produce .hex files for flashing in `../bld/<app>/fabric/<app>.hex`

### Flashing

To flash an app to the EVK, make sure your BOOT switches on the board are set to `101` and run the following command from an app's `fabric` folder where the .hex is contained:

```
eff-flash <app>.hex sram
```

If you'd like to flash the app to non-volatile memory, please change `sram` to `mram`, flash, power off your board, set BOOT switches to 010, then power on to run the app.

### Monitoring

To inspect the app output please run the following command from the terminal:

```
minicom -b 115200 -D /dev/ttyACM2
```
