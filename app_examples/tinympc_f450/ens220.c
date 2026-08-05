/**
 ******************************************************************************
 *
 * @file ens220.c
 *
 * @brief ScioSense ENS220 Barometric Pressure and Temperature Sensor driver
 *
 * The ENS220 uses standard 8-bit-register I2C (eff_i2c_read / eff_i2c_write),
 * the same interface as BMM350.
 *
 * Measurement output registers (read as one 5-byte burst from 0x17):
 *   [PRESS_LSB][PRESS_MID][PRESS_MSB][TEMP_LSB][TEMP_MSB]
 *
 * Conversion (integer-only):
 *   P[hPa] = raw_p / 6400        (frac = (raw_p % 6400) * 100 / 6400)
 *   T[°C]  = (raw_t − 34963) / 128  (34963 = round(273.15 × 128))
 *
 ******************************************************************************
 */

/*
 * INCLUDE FILES
 ******************************************************************************
 */
#include "ens220.h"
#include <eff.h>
#include <eff/drivers/i2c.h>
#include <stdio.h>
#include <stdint.h>

/*
 * EXTERNAL DECLARATIONS
 ******************************************************************************
 */

/* I2C peripheral shared by all hat sensors — Bank 2 (DIGIO024=SCL, DIGIO025=SDA) */
extern eff_i2c_t *I2C_2_1;

/*
 * PUBLIC FUNCTIONS
 ******************************************************************************
 */

int8_t ens220_init(void)
{
    uint8_t val;

    eff_pinmux_set(PINMUX_2, PINMUX_SPI_I2C1);
    eff_i2c_init(I2C_2_1, I2C_SPEED_100K);
    eff_gpio_pull_set(GPIO_2, GPIO_PIN_4 | GPIO_PIN_5, EFF_GPIO_PULL_NONE);

    /* Soft reset */
    val = ENS220_MODE_RESET;
    eff_i2c_write(I2C_2_1, ENS220_I2C_ADDR, ENS220_REG_MODE_CFG, &val, 1);
    for (volatile uint32_t i = 0; i < ENS220_RESET_CYCLES; i++) {}

    /* One-shot standby: auto-idle after each measurement */
    val = ENS220_STBY_ONESHOT;
    eff_i2c_write(I2C_2_1, ENS220_I2C_ADDR, ENS220_REG_STBY_CFG, &val, 1);

    printf("ens220_init: OK\n");
    return 0;
}

int8_t ens220_read_id(void)
{
    volatile uint8_t buf[2] = {0};
    uint16_t id;

    eff_i2c_read(I2C_2_1, ENS220_I2C_ADDR, ENS220_REG_PART_ID, (uint8_t *)buf, 2);
    id = (uint16_t)buf[0] | ((uint16_t)buf[1] << 8);   /* little-endian */

    if (id != ENS220_PART_ID_EXPECTED) {
        printf("ens220_id: 0x%04X — unexpected\n", (unsigned int)id);
        return -1;
    }

    printf("ens220_id: 0x%04X — ENS220 confirmed\n", (unsigned int)id);
    return 0;
}

int8_t ens220_reset(void)
{
    uint8_t val = ENS220_MODE_RESET;

    eff_i2c_write(I2C_2_1, ENS220_I2C_ADDR, ENS220_REG_MODE_CFG, &val, 1);
    for (volatile uint32_t i = 0; i < ENS220_RESET_CYCLES; i++) {}

    printf("ens220_reset: OK\n");
    return 0;
}

int8_t ens220_read_raw(uint32_t *raw_p, uint16_t *raw_t)
{
    uint8_t val;
    volatile uint8_t status = 0;
    volatile uint8_t buf[5] = {0};

    val = ENS220_MODE_MEASURE;
    eff_i2c_write(I2C_2_1, ENS220_I2C_ADDR, ENS220_REG_MODE_CFG, &val, 1);

    for (int tries = 0; tries < 50; tries++) {
        for (volatile uint32_t i = 0; i < ENS220_POLL_CYCLES; i++) {}
        eff_i2c_read(I2C_2_1, ENS220_I2C_ADDR, ENS220_REG_DATA_STAT, (uint8_t *)&status, 1);
        if ((status & 0x03u) == 0x03u)
            break;
    }

    if ((status & 0x03u) != 0x03u)
        return -1;

    eff_i2c_read(I2C_2_1, ENS220_I2C_ADDR, ENS220_REG_PRESS_OUT, (uint8_t *)buf, 5);

    *raw_p = (uint32_t)buf[0] | ((uint32_t)buf[1] << 8) | ((uint32_t)buf[2] << 16);
    *raw_t = (uint16_t)buf[3] | ((uint16_t)buf[4] << 8);
    return 0;
}

int8_t ens220_measure(void)
{
    uint32_t raw_p;
    uint16_t raw_t;

    if (ens220_read_raw(&raw_p, &raw_t) != 0) {
        printf("ens220_measure: timeout\n");
        return -1;
    }

    uint32_t p_hpa  = raw_p / 6400;
    uint32_t p_frac = (raw_p % 6400) * 100 / 6400;
    int32_t  t_off  = (int32_t)raw_t - 34963;
    int32_t  t_int  = t_off / 128;
    int32_t  t_abs  = t_off < 0 ? -t_off : t_off;
    int32_t  t_frac = t_abs % 128 * 10 / 128;

    printf("ens220: raw_p=%u raw_t=%u P=%lu.%02lu hPa",
           (unsigned)raw_p, (unsigned)raw_t,
           (unsigned long)p_hpa, (unsigned long)p_frac);
    if (t_int == 0 && t_off < 0)
        printf(" T=-0.%ld C\n", (long)t_frac);
    else
        printf(" T=%ld.%ld C\n", (long)t_int, (long)t_frac);

    return 0;
}
