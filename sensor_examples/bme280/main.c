#include <eff.h>
#include <stdio.h>
#include "bme280.h"

/* E1x Configuration */
#define PERIPH_MODE             PINMUX_I2C0_I2C1
#define PINMUX_PORT             PINMUX_0
#define I2C_DEVICE              I2C_0_0
#define I2C_SPEED               I2C_SPEED_100K

// BME280 I2C address (0x76 or 0x77, depending on SDO pin)
#define BME280_I2C_ADDR 0x77

// BME280 API requires these I2C interface functions
int8_t user_i2c_read(uint8_t reg_addr, uint8_t *reg_data, uint32_t len, void *intf_ptr)
{
    // printf("About to perform eff_i2c_read...\r\n");
    // printf("reg_addr: 0x%02X\r\n", reg_addr);
    // printf("len: %d\r\n", len);
    eff_i2c_t* i2c = (eff_i2c_t*)intf_ptr;
    int8_t result = eff_i2c_read(i2c, BME280_I2C_ADDR, reg_addr, reg_data, len);
    // printf("eff_i2c_read result: %d\r\n", result);
    // printf("reg_data: 0x%02X\r\n", *reg_data);
    return result;
}

int8_t user_i2c_write(uint8_t reg_addr, const uint8_t *reg_data, uint32_t len, void *intf_ptr)
{
    eff_i2c_t* i2c = (eff_i2c_t*)intf_ptr;
    return eff_i2c_write(i2c, BME280_I2C_ADDR, reg_addr, (uint8_t*)reg_data, len);
}

void user_delay_us(uint32_t period, void *intf_ptr)
{
    (void)intf_ptr;
    sleep_us(period);
}

int main(void)
{
    printf("Starting BME280 sensor application...\r\n");
    struct bme280_dev dev;
    int8_t rslt;
    uint8_t chip_id = 0;
    struct bme280_data comp_data; // Compensated data from the sensor

    // Configure pinmux for I2C
    eff_pinmux_set(PINMUX_PORT, PERIPH_MODE);
    
    // Initialize I2C
    eff_i2c_init(I2C_DEVICE, I2C_SPEED);

    // Initialize BME280 device structure
    // This section maps the BME280 API to the Efficient SDK I2C functions
    dev.intf = BME280_I2C_INTF;
    dev.intf_ptr = I2C_DEVICE;
    dev.read = user_i2c_read;
    dev.write = user_i2c_write;
    dev.delay_us = user_delay_us;

    // Initialize BME280
    rslt = bme280_init(&dev);
    if (rslt != BME280_OK) {
        printf("BME280 init failed: %d\r\n", rslt);
        return -1;
    }

    // Read chip ID (stored in dev.chip_id after init, but read explicitly)
    rslt = bme280_get_regs(BME280_REG_CHIP_ID, &chip_id, 1, &dev);
    if (rslt != BME280_OK) {
        printf("Failed to read chip ID: %d\r\n", rslt);
        return -1;
    }
    printf("BME280 Chip ID: 0x%02X\r\n", chip_id);

    // Configure sensor settings
    // Added oversampling may be necessary for better accuracy
    struct bme280_settings settings;
    settings.osr_h = BME280_OVERSAMPLING_1X;  // Humidity oversampling
    settings.osr_p = BME280_OVERSAMPLING_1X;  // Pressure oversampling, extra oversampling adds accuracy
    settings.osr_t = BME280_OVERSAMPLING_1X;  // Temperature oversampling
    // settings.filter = BME280_FILTER_COEFF_16;  // Filter coefficient 16 for noise reduction
    settings.filter = BME280_FILTER_COEFF_OFF;
    settings.standby_time = BME280_STANDBY_TIME_0_5_MS;
    
    rslt = bme280_set_sensor_settings(BME280_SEL_OSR_PRESS | BME280_SEL_OSR_TEMP | BME280_SEL_OSR_HUM | BME280_SEL_FILTER | BME280_SEL_STANDBY, &settings, &dev);
    if (rslt != BME280_OK) {
        printf("Failed to set sensor settings: %d\r\n", rslt);
        return -1;
    }

    while (1) {
        // Set sensor to forced mode to trigger a measurement
        // A single measurement is performed, the sensor goes back to sleep mode,
        // and the measurement results can be obtained from data registers
        rslt = bme280_set_sensor_mode(BME280_POWERMODE_FORCED, &dev);
        if (rslt != BME280_OK) {
            printf("Failed to set sensor mode: %d\r\n", rslt);
            return -1;
        }

        // Wait for measurement to complete (calculates delay based on settings)
        uint32_t req_delay;
        bme280_cal_meas_delay(&req_delay, &settings);
        dev.delay_us(req_delay, dev.intf_ptr);

        // Read temperature, pressure, and humidity
        rslt = bme280_get_sensor_data(BME280_ALL, &comp_data, &dev);
        if (rslt != BME280_OK) {
            printf("Failed to read sensor data: %d\r\n", rslt);
            return -1;
        }

        // Temperature
        // Temperature is in 0.01°C units (e.g., 2500 = 25.00°C)
        int32_t temp = comp_data.temperature;
        int32_t temp_whole = temp / 100;
        int32_t temp_frac = (temp < 0 ? -temp : temp) % 100;
        printf("Temperature: %d.%02d C\r\n", temp_whole, temp_frac);
    
        // Pressure
        // In 32-bit mode, pressure is directly in Pa
        int32_t pressure_pa = (int32_t)comp_data.pressure;
        printf("Pressure: %d Pa\r\n", pressure_pa);
    
        // Humidity
        // The integer mode returns humidity in units where 102400 = 100.00%
        // This means 1% = 1024 units, so we divide by 1024 to get percentage
        uint32_t hum = comp_data.humidity;
        uint32_t hum_whole = hum / 1024; // Convert from 1/1024% units to percentage
        uint32_t hum_frac = ((hum % 1024) * 100) / 1024;
        printf("Humidity: %lu.%02lu %%\r\n", (unsigned long)hum_whole, (unsigned long)hum_frac);

        printf("\r\n");

        sleep_ms(1000);
    }

    return 0;
}

