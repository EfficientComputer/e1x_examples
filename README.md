# E1x Examples

Welcome to the E1x Examples repository! This collection of example applications is designed to help you get up and running quickly with the E1x processor and the E1x Evaluation Kit (EVK). Whether you're exploring E1x for the first time or building production applications, these examples provide a solid starting point and reference implementations that demonstrate best practices for developing on Efficient hardware. Use them as templates for your own projects, or as building blocks to accelerate your development.

If you haven't already, please make sure you've set up your board and development environment using our [Evaluation Kit Setup Instructions](https://docs.efficient.computer/evaluation-kit).

---

## App Examples

Optimized algorithm implementations demonstrating how to leverage E1x's architecture for real-world DSP, image processing, communications, and compute workloads. Each app includes complete source code with E1x-optimized kernels and correctness validation. This project supports both our physical EVK and our Cloud EVK.

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
- **tinympc** - TinyMPC Quadrotor Model-Predictive Control (Robotics / Control)
- **tinympc_lpv** - Cache-Scheduled (LPV) Fixed-Wing MPC (Robotics / Control)
- **mekf** - MEKF IMU Attitude Estimator (Sensor Fusion / Robotics)
- **quickstart** - Simple "Hello World" example
---

## Sensor Examples

Example apps demonstrating how to interface various sensor modules with the E1x processor. These examples show sensor initialization, configuration, and data acquisition using the E1x SDK's peripheral drivers. Please note that this project does not support our Cloud EVK.

- **adxl372** - ADXL372 High-G Accelerometer
- **bme280** - BME280 Environmental Sensor (Temperature, Humidity, Pressure)
- **arducam_ov2640** - Arducam Mini 2MP Plus Camera Module
