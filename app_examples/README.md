# E1x App Examples

## Overview

This project provides example apps demonstrating optimized algorithms running on the E1x processor using the E1x Evaluation Kit (EVK). These examples serve as a reference for developers looking to understand how to leverage E1x's architecture for real-world DSP, image processing, communications, and compute workloads. Each app includes complete source code showing EFF-optimized implementations along with correctness validation.

---

## Apps

* **biquad_filter** - Optimized Biquad IIR Filter (DSP)
* **conv5x5** - Optimized 5×5 Convolution Kernel (ML/Vision)
* **dmm** - Dense Matrix Multiply (Linear Algebra)
* **fft** - Fast Fourier Transform (DSP)
* **jpeg** - JPEG Image Encoder (Multimedia)
* **ldpc** - LDPC Encoder/Decoder (Communications)
* **mlperftiny** - MLPerfTiny neural network benchmark suite
* **quickstart** - Simple "Hello World" example

## Getting Started

If you haven't already, please make sure you've set up your board and development environment using our [Evaluation Kit Setup Instructions](https://docs.efficient.computer/evaluation-kit). This project supports both our physical EVK and our Cloud EVK.

Additional details are contained within each app's folder.

The **mlperftiny** folder contains five MLPerfTiny benchmark apps:
anomaly detection, image classification, keyword spotting, streaming wake word,
and visual wake words. Each benchmark has its own README with model, input,
target, and validation details.

### App Configuration

Each app is contained in its own directory with a `CMakeLists.txt` for building. Refer to the individual app README and source files for detailed usage and configuration options.

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
