# BME280 - Environmental Sensor

## Overview

This sensor can be purchased at [Digikey](https://www.digikey.com/short/dz5r4cqm).

The BME280 example demonstrates how to interface with Bosch's BME280 environmental sensor over I2C. This sensor captures three key environmental parameters: **temperature**, **barometric pressure**, and **relative humidity**. It leverages Bosch's Sensor API as it contains crucial configuration options and capabilities. The example shows how to configure the sensor's oversampling and filter settings, trigger measurements in forced mode, and read the compensated sensor data. This makes it ideal for weather monitoring stations, indoor climate control systems, and environmental data logging applications.

| Measurement | Unit |
|-------------|------|
| Temperature | °C |
| Pressure | Pa |
| Humidity | % RH |

---

## Connecting Your Sensor

Connect your sensor module to the EVK using the pin mappings shown below. Please note that Bank 0 is used for the BME280 and Bank 2 is used for the ADXL372. This was done to show the use of different peripheral modes on different Banks and the pin mappings associated with them. For more information peripheral modes and pin mappings please see Sections 7.3 and 7.4 of the [E1x Datasheet](https://docs.efficient.computer/assets/files/E1M3M4M-Datasheet-v0.10-505a7892ef5fb5e52b9c52952fb9e1c2.pdf).

Note: The pin mappings are shown in first two columns. The second and third columns have been added for additional info/context.

### BME280 (I2C — Bank 0, Mode 9, Datasheet 7.4)

| BME280 | EVK | DIGIO    | Peripheral   |
|--------|-----|----------|--------------|
| VIN    | 3V3 | —        | —            |
| GND    | GND | —        | —            |
| SCK    | 000 | DIGIO000 | I2C0_SCL [0] |
| SDI    | 001 | DIGIO001 | I2C0_SDA [1] |