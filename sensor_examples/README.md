# E1x Sensor Examples

## Overview

This repository provides example apps demonstrating how to interface various sensor modules with the E1x processor using the E1x Evaluation Kit (EVK). These examples serve as a reference for developers looking to integrate real-world sensing capabilities into their E1x-based applications. Each app includes complete source code showing sensor initialization, configuration, and data acquisition using the E1x SDK's peripheral drivers.

---

## Apps

### BME280 — Environmental Sensor

Can be purchased at [Digikey](https://www.digikey.com/short/dz5r4cqm).

The BME280 example demonstrates how to interface with Bosch's BME280 environmental sensor over I2C. This sensor captures three key environmental parameters: **temperature**, **barometric pressure**, and **relative humidity**. It leverages Bosch's Sensor API as it contains crucial configuration options and capabilities. The example shows how to configure the sensor's oversampling and filter settings, trigger measurements in forced mode, and read the compensated sensor data. This makes it ideal for weather monitoring stations, indoor climate control systems, and environmental data logging applications.

| Measurement | Unit |
|-------------|------|
| Temperature | °C |
| Pressure | Pa |
| Humidity | % RH |

---

### ADXL372 — High-G Accelerometer

Can be purchased at [Digikey](https://www.digikey.com/short/wp5q9349).

The ADXL372 example demonstrates how to interface with Analog Devices' ADXL372 3-axis accelerometer over SPI. While an API exists for this sensor, it was not used given the relative simplicity of accessing readings and the desire to preserve the simplicity of the example. This sensor is designed for high-g shock and vibration detection, capable of measuring acceleration up to ±200g with 12-bit resolution. The example covers device initialization, register configuration, and continuous reading of acceleration data on the **X**, **Y**, and **Z** axes. This sensor is well-suited for impact detection, structural health monitoring, and motion sensing in demanding environments.

| Measurement | Scale | Resolution |
|-------------|-------|------------|
| X-axis acceleration | 100 mg/LSB | 12-bit |
| Y-axis acceleration | 100 mg/LSB | 12-bit |
| Z-axis acceleration | 100 mg/LSB | 12-bit |

---

## Getting Started

If you haven't already, please make sure you've set up your board and development environment using our [Evaluation Kit Setup Instructions](https://docs.efficient.computer/evaluation-kit).

### Connecting Your Sensor

Connect your sensor module to the EVK using the pin mappings shown below. Please note that Bank 0 is used for the BME280 and Bank 2 is used for the ADXL372. This was done to show the use of different peripheral modes on different Banks and the pin mappings associated with them. For more information peripheral modes and pin mappings please see Sections 7.3 and 7.4 of the [E1x Datasheet](https://docs.efficient.computer/assets/files/E1M3M4M-Datasheet-v0.10-505a7892ef5fb5e52b9c52952fb9e1c2.pdf).

Note: The pin mappings are shown in first two columns. The second and third columns have been added for additional info/context.

#### BME280 (I2C — Bank 0, Mode 9, Datasheet 7.4)

| BME280 | EVK | DIGIO    | Peripheral   |
|--------|-----|----------|--------------|
| VIN    | 3V3 | —        | —            |
| GND    | GND | —        | —            |
| SCK    | 000 | DIGIO000 | I2C0_SCL [0] |
| SDI    | 001 | DIGIO001 | I2C0_SDA [1] |

#### ADXL372 (SPI — Bank 2, Mode 3, Datasheet 7.4)

| ADXL372 | EVK | DIGIO    | Peripheral   |
|---------|-----|----------|--------------|
| VS      | 1V8 | —        | —            |
| VIO     | 1V8 | —        | —            |
| GND     | GND | —        | —            |
| SCLK    | 020 | DIGIO020 | SPI_CLK [0]  |
| CS      | 021 | DIGIO021 | SPI_CS [1]   |
| MOSI    | 022 | DIGIO022 | SPI_MOSI [2] |
| MISO    | 023 | DIGIO023 | SPI_MISO [3] |

### App/Sensor Configuration

Each app is contained in its own directory with a `CMakeLists.txt` for building. Refer to the individual app source files for detailed usage and configuration options for each sensor.

### Building

To build all apps, execute the following commands from the top level of this folder:

```
mkdir build
cd build
cmake -G Ninja .. -DEFF_STDIO_PORT=3
ninja
```

This will produce .hex files for flashing in `e1x_sensor_examples/build/<app>/scalar/<app>.hex`

### Flashing

To flash an app to the EVK, make sure your BOOT switches on the board are set to 101 and run the following command from an app's scalar folder where the .hex is contained:

```
eff-flash <app>.hex sram
```

If you'd like to flash the app to non-volatile memory, please change `sram` to `mram`, flash, power off your board, set BOOT switches to 010, then power on to run the app.

### Monitoring

To inspect the sensor output please run the following command from the terminal:

```
minicom -b 115200 -D /dev/ttyACM2
```
