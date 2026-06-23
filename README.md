# E1x Examples

Welcome to the E1x Examples repository! This collection of example applications is designed to help you get up and running quickly with the E1x processor and the E1x Evaluation Kit (EVK). Whether you're exploring E1x for the first time or building production applications, these examples provide a solid starting point and reference implementations that demonstrate best practices for developing on EFF hardware. Use them as templates for your own projects, or as building blocks to accelerate your development.

If you haven't already, please make sure you've set up your board and development environment using our [Evaluation Kit Setup Instructions](https://docs.efficient.computer/evaluation-kit).

---

## App Examples

Optimized algorithm implementations demonstrating how to leverage E1x's architecture for real-world DSP, image processing, communications, and compute workloads. Each app includes complete source code with EFF-optimized kernels and correctness validation. This project supports both our physical EVK and our Cloud EVK.

* **biquad_filter** - Optimized Biquad IIR Filter (DSP)
* **conv5x5** - Optimized 5×5 Convolution Kernel (ML/Vision)
* **dmm** - Dense Matrix Multiply (Linear Algebra)
* **fft** - Fast Fourier Transform (DSP)
* **jpeg** - JPEG Image Encoder (Multimedia)
* **ldpc** - LDPC Encoder/Decoder (Communications)
* **tinympc** - TinyMPC Quadrotor Model-Predictive Control (Robotics / Control)
* **tinympc_lpv** - Cache-Scheduled (LPV) Fixed-Wing MPC (Robotics / Control)
* **mekf** - MEKF IMU Attitude Estimator (Sensor Fusion / Robotics)
* **quickstart** - Simple "Hello World" example

---

## Sensor Examples

Example apps demonstrating how to interface various sensor modules with the E1x processor. These examples show sensor initialization, configuration, and data acquisition using the E1x SDK's peripheral drivers. Please note that this project does not support our Cloud EVK.

* **adxl372** - ADXL372 High-G Accelerometer
* **bme280** - BME280 Environmental Sensor (Temperature, Humidity, Pressure)
