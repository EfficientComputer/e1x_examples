/**
 ******************************************************************************
 *
 * @file bmi323.c
 *
 * @brief Bosch BMI323 6DoF IMU driver
 *
 * Implements I2C register access for the BMI323 accelerometer + gyroscope
 * on the Efficient Computer E1x rev3 EVK.
 *
 * Hardware quirk: BMI323 inserts 2 dummy bytes before payload on every
 * I2C read (datasheet section 4). All reads request 4 bytes and print
 * the full buffer to determine correct byte ordering.
 *
 ******************************************************************************
 */

/*
 * INCLUDE FILES
 ******************************************************************************
 */
#include "bmi323.h"
#include <eff.h>
#include <eff/drivers/i2c.h>
#include <stdio.h>
#include <stdint.h>

/*
 * DEFINES
 ******************************************************************************
 */

/* I2C peripheral shared by all hat sensors — Bank 2 (DIGIO024=SCL, DIGIO025=SDA) */
extern eff_i2c_t *I2C_2_1;

/*
 * Spin-loop delay for POR/suspend idle time.
 * BMI323 requires >= 450 µs between accesses after reset.
 * At 50 MHz CPU: 25000 cycles ≈ 500 µs.
 */
#define BMI323_POR_DELAY_CYCLES  25000u

/*
 * REGISTER ACCESS
 ******************************************************************************
 */

int8_t bmi323_reg_read(uint8_t reg, uint16_t *val)
{
    volatile uint8_t buf[4] = {0};
    int8_t rc = eff_i2c_read(I2C_2_1, BMI323_I2C_ADDR, reg, buf, 4);
    if (rc != 0) {
        printf("bmi323: read 0x%02x failed\n", reg);
        return -1;
    }
    /* buf[0..1] are dummy bytes; payload is buf[2..3] (little-endian) */
    *val = (uint16_t)buf[2] | ((uint16_t)buf[3] << 8);
    return 0;
}

int8_t bmi323_reg_write(uint8_t reg, uint16_t val)
{
    uint8_t buf[2] = { (uint8_t)(val & 0xFF), (uint8_t)(val >> 8) };
    int8_t rc = eff_i2c_write(I2C_2_1, BMI323_I2C_ADDR, reg, buf, 2);
    if (rc != 0) {
        printf("bmi323: write 0x%02x failed\n", reg);
        return -1;
    }
    return 0;
}

/*
 * PUBLIC FUNCTIONS
 ******************************************************************************
 */

int8_t bmi323_init(void)
{
    eff_pinmux_set(PINMUX_2, PINMUX_SPI_I2C1);

    if (eff_i2c_init(I2C_2_1, I2C_SPEED_100K) != 0) {
        printf("bmi323_init: I2C init failed\n");
        return -1;
    }
    eff_gpio_pull_set(GPIO_2, GPIO_PIN_4 | GPIO_PIN_5, EFF_GPIO_PULL_NONE);

    printf("bmi323_init: OK\n");
    return 0;
}

int8_t bmi323_check_id(void)
{
    uint16_t val = 0;
    if (bmi323_reg_read(BMI323_REG_CHIP_ID, &val) != 0)
        return -1;

    if ((val & 0xFF) != BMI323_CHIP_ID_EXPECTED) {
        printf("bmi323_check_id: unexpected CHIP_ID 0x%02x", (unsigned int)(val & 0xFF));
        printf(" (expected 0x%02x)\n", (unsigned int)BMI323_CHIP_ID_EXPECTED);
        return -1;
    }

    printf("bmi323_check_id: CHIP_ID = 0x%02x — BMI323 confirmed\n",
           (unsigned int)(val & 0xFF));
    return 0;
}

int8_t bmi323_reset(void)
{
    if (bmi323_reg_write(BMI323_REG_CMD, BMI323_CMD_SOFT_RESET) != 0)
        return -1;

    sleep_ms(5);  /* POR sequence: datasheet min 450 µs, use 5 ms for margin */

    printf("bmi323_reset: done\n");
    return 0;
}

int8_t bmi323_configure(void)
{
    uint16_t err = 0, acc_rb = 0, gyr_rb = 0, status = 0;

    /* Read ERR_REG before configuring — non-zero means sensor entered a fault state */
    bmi323_reg_read(BMI323_REG_ERR_REG, &err);
    printf("bmi323_configure: ERR_REG = 0x%04x\n", (unsigned int)err);

    /* Accelerometer: normal mode, 100 Hz, ±2g */
    if (bmi323_reg_write(BMI323_REG_ACC_CONF, BMI323_ACC_CONF_NORMAL) != 0)
        return -1;
    sleep_ms(1);

    /* Gyroscope: normal mode, 100 Hz, ±2000°/s */
    if (bmi323_reg_write(BMI323_REG_GYR_CONF, BMI323_GYR_CONF_NORMAL) != 0)
        return -1;
    sleep_ms(1);

    /* Read back ACC_CONF and GYR_CONF to verify the writes landed */
    bmi323_reg_read(BMI323_REG_ACC_CONF, &acc_rb);
    bmi323_reg_read(BMI323_REG_GYR_CONF, &gyr_rb);
    printf("bmi323_configure: ACC_CONF = 0x%04x (wrote 0x%04x) — %s\n",
           (unsigned int)acc_rb, (unsigned int)BMI323_ACC_CONF_NORMAL,
           acc_rb == BMI323_ACC_CONF_NORMAL ? "OK" : "MISMATCH");
    printf("bmi323_configure: GYR_CONF = 0x%04x (wrote 0x%04x) — %s\n",
           (unsigned int)gyr_rb, (unsigned int)BMI323_GYR_CONF_NORMAL,
           gyr_rb == BMI323_GYR_CONF_NORMAL ? "OK" : "MISMATCH");

    if (acc_rb != BMI323_ACC_CONF_NORMAL || gyr_rb != BMI323_GYR_CONF_NORMAL) {
        printf("bmi323_configure: config write failed — sensor not ready for configuration?\n");
        return -1;
    }

    /* Wait for first sample at 100 Hz; use 100 ms for margin */
    sleep_ms(100);

    /* STATUS[7]=drdy_acc, STATUS[6]=drdy_gyr, STATUS[4]=cmd_rdy */
    bmi323_reg_read(BMI323_REG_STATUS, &status);
    printf("bmi323_configure: STATUS = 0x%02x  (drdy_acc=%d drdy_gyr=%d)\n",
           (unsigned int)(status & 0xFF),
           (int)((status >> 7) & 1),
           (int)((status >> 6) & 1));

    if (!((status >> 7) & 1) || !((status >> 6) & 1)) {
        printf("bmi323_configure: sensor did not enter normal mode\n");
        return -1;
    }

    printf("bmi323_configure: normal mode, 100 Hz, accel ±2g, gyro ±2000°/s\n");
    return 0;
}

/*
 * BMI323 data registers 0x03-0x08 must be read in a single burst transaction.
 * Individual per-register reads return 0. Both read_accel and read_gyro burst
 * all 7 data registers (accel + gyro + temp) in one 16-byte read (2 platform
 * dummies + 14 data bytes) and extract the relevant axes.
 *
 * Burst layout:
 *   buf[0..1]   = platform dummies
 *   buf[2..3]   = ACC_X    buf[4..5]  = ACC_Y    buf[6..7]  = ACC_Z
 *   buf[8..9]   = GYR_X    buf[10..11]= GYR_Y    buf[12..13]= GYR_Z
 *   buf[14..15] = TEMP
 */
int8_t bmi323_read_accel(int16_t *x, int16_t *y, int16_t *z)
{
    uint16_t status = 0;
    volatile uint8_t buf[16] = {0};

    /* STATUS read is required to latch the data registers before reading. */
    bmi323_reg_read(BMI323_REG_STATUS, &status);

    if (eff_i2c_read(I2C_2_1, BMI323_I2C_ADDR,
                     BMI323_REG_ACC_DATA_X, (uint8_t *)buf, 16) != 0)
        return -1;

    *x = (int16_t)((uint16_t)buf[2]  | ((uint16_t)buf[3]  << 8));
    *y = (int16_t)((uint16_t)buf[4]  | ((uint16_t)buf[5]  << 8));
    *z = (int16_t)((uint16_t)buf[6]  | ((uint16_t)buf[7]  << 8));
    return 0;
}

int8_t bmi323_read_gyro(int16_t *x, int16_t *y, int16_t *z)
{
    uint16_t status = 0;
    volatile uint8_t buf[16] = {0};

    bmi323_reg_read(BMI323_REG_STATUS, &status);

    if (eff_i2c_read(I2C_2_1, BMI323_I2C_ADDR,
                     BMI323_REG_ACC_DATA_X, (uint8_t *)buf, 16) != 0)
        return -1;

    *x = (int16_t)((uint16_t)buf[8]  | ((uint16_t)buf[9]  << 8));
    *y = (int16_t)((uint16_t)buf[10] | ((uint16_t)buf[11] << 8));
    *z = (int16_t)((uint16_t)buf[12] | ((uint16_t)buf[13] << 8));
    return 0;
}

int8_t bmi323_read_temp(int16_t *temp)
{
    uint16_t raw = 0;

    if (bmi323_reg_read(BMI323_REG_TEMP_DATA, &raw) != 0)
        return -1;

    *temp = (int16_t)raw;
    return 0;
}
