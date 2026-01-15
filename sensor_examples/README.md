# E1x Sensor Examples

## Overview

This project provides example apps demonstrating how to interface various sensor modules with the E1x processor using the E1x Evaluation Kit (EVK). These examples serve as a reference for developers looking to integrate real-world sensing capabilities into their E1x-based applications. Each app includes complete source code showing sensor initialization, configuration, and data acquisition using the E1x SDK's peripheral drivers.

---

## Apps

* ADXL372 - High-G Accelerometer
* BME280 - Environmental Sensor

## Getting Started

If you haven't already, please make sure you've set up your board and development environment using our [Evaluation Kit Setup Instructions](https://docs.efficient.computer/evaluation-kit). Please note that this this project does not support our Cloud EVK.

Additional details as well as connection instructions are contained within each app's folder.


### App/Sensor Configuration

Each app is contained in its own directory with a `CMakeLists.txt` for building. Refer to the individual app source files for detailed usage and configuration options for each sensor.

### Building

To build all apps, execute the following commands from the top level of this folder:

```
mkdir bld
cd bld
cmake -G Ninja .. -DEFF_STDIO_PORT=3
ninja
```

This will produce .hex files for flashing in `../bld/<app>/scalar/<app>.hex`
Please note that when you reach the point of running AI workloads on E1x, you'll build for `fabric` instead of `scalar`. That setting is adjusted in each app's `CMakeLists.txt`.

### Flashing

To flash an app to the EVK, make sure your BOOT switches on the board are set to `101` and run the following command from an app's scalar folder where the .hex is contained:

```
eff-flash <app>.hex sram
```

If you'd like to flash the app to non-volatile memory, please change `sram` to `mram`, flash, power off your board, set BOOT switches to 010, then power on to run the app.

### Monitoring

To inspect the sensor output please run the following command from the terminal:

```
minicom -b 115200 -D /dev/ttyACM2
```
