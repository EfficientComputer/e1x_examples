# ADXL372 — High-G Accelerometer

## Overview

This sensor can be purchased at [Digikey](https://www.digikey.com/short/wp5q9349).

The ADXL372 example demonstrates how to interface with Analog Devices' ADXL372 3-axis accelerometer over SPI. While an API exists for this sensor, it was not used given the relative simplicity of accessing readings and the desire to preserve the simplicity of the example. This sensor is designed for high-g shock and vibration detection, capable of measuring acceleration up to ±200g with 12-bit resolution. The example covers device initialization, register configuration, and continuous reading of acceleration data on the **X**, **Y**, and **Z** axes. This sensor is well-suited for impact detection, structural health monitoring, and motion sensing in demanding environments.

| Measurement | Scale | Resolution |
|-------------|-------|------------|
| X-axis acceleration | 100 mg/LSB | 12-bit |
| Y-axis acceleration | 100 mg/LSB | 12-bit |
| Z-axis acceleration | 100 mg/LSB | 12-bit |

---

## Connecting Your Sensor

Connect your sensor module to the EVK using the pin mappings shown below. Please note that Bank 0 is used for the BME280 and Bank 2 is used for the ADXL372. This was done to show the use of different peripheral modes on different Banks and the pin mappings associated with them. For more information peripheral modes and pin mappings please see Sections 7.3 and 7.4 of the [E1x Datasheet](https://docs.efficient.computer/assets/files/E1M3M4M-Datasheet-v0.10-505a7892ef5fb5e52b9c52952fb9e1c2.pdf).

Note: The pin mappings are shown in first two columns. The second and third columns have been added for additional info/context.

### ADXL372 (SPI — Bank 2, Mode 3, Datasheet 7.4)

| ADXL372 | EVK | DIGIO    | Peripheral   |
|---------|-----|----------|--------------|
| VS      | 1V8 | —        | —            |
| VIO     | 1V8 | —        | —            |
| GND     | GND | —        | —            |
| SCLK    | 020 | DIGIO020 | SPI_CLK [0]  |
| CS      | 021 | DIGIO021 | SPI_CS [1]   |
| MOSI    | 022 | DIGIO022 | SPI_MOSI [2] |
| MISO    | 023 | DIGIO023 | SPI_MISO [3] |