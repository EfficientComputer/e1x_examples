/**
 ******************************************************************************
 *
 * @file bmi323.h
 *
 * @brief Bosch BMI323 6DoF IMU driver — public API
 *
 * BMI323 accelerometer + gyroscope
 * Interface: I2C (8-bit register addresses, 16-bit register values)
 * I2C address: 0x68
 * Host: Efficient Computer E1x rev3 EVK
 *
 * NOTE: BMI323 inserts 2 dummy bytes before payload on every I2C read.
 * The driver must account for these dummy bytes when parsing the read buffer.
 *
 ******************************************************************************
 */

#ifndef _BMI323_H_
#define _BMI323_H_

/*
 * INCLUDE FILES
 ******************************************************************************
 */
#include <eff.h>
#include <eff/drivers/i2c.h>
#include <stdint.h>

/*
 * MACROS — I2C
 ******************************************************************************
 */
#define BMI323_I2C_ADDR          0x68u

/*
 * MACROS — register addresses
 ******************************************************************************
 */
#define BMI323_REG_CHIP_ID           0x00u
#define BMI323_REG_ERR_REG           0x01u
#define BMI323_REG_STATUS            0x02u
#define BMI323_REG_ACC_DATA_X        0x03u
#define BMI323_REG_ACC_DATA_Y        0x04u
#define BMI323_REG_ACC_DATA_Z        0x05u
#define BMI323_REG_GYR_DATA_X        0x06u
#define BMI323_REG_GYR_DATA_Y        0x07u
#define BMI323_REG_GYR_DATA_Z        0x08u
#define BMI323_REG_TEMP_DATA         0x09u
#define BMI323_REG_INT_STATUS_INT1   0x0Du
#define BMI323_REG_INT_STATUS_INT2   0x0Eu
#define BMI323_REG_ACC_CONF          0x20u
#define BMI323_REG_GYR_CONF          0x21u
#define BMI323_REG_CMD               0x7Eu

/*
 * MACROS — configuration values
 ******************************************************************************
 */

/* ACC_CONF: acc_mode=0b011(low-power duty-cycle), odr=0b0001(0.78 Hz),
 *           range=0b010(±8g), avg=0, bw=0 */
#define BMI323_ACC_CONF_LOWPWR   0x3021u

/* GYR_CONF: gyr_mode=0b011(low-power duty-cycle), odr=0b0001(0.78 Hz),
 *           range=0b100(±2000°/s), avg=0, bw=0 */
#define BMI323_GYR_CONF_LOWPWR   0x3041u

/* ACC_CONF: acc_mode=0b100(normal), odr=0b1000(100 Hz), range=0b000(±2g), avg=0, bw=0 */
#define BMI323_ACC_CONF_NORMAL   0x4008u

/* GYR_CONF: gyr_mode=0b100(normal), odr=0b1000(100 Hz), range=0b100(±2000°/s), avg=0, bw=0 */
#define BMI323_GYR_CONF_NORMAL   0x4048u

/* Sensitivity at ±2g range: 2^15 / 2 = 16384 counts per g */
#define BMI323_ACC_COUNTS_PER_G  16384

#define BMI323_CHIP_ID_EXPECTED  0x43u
#define BMI323_CMD_SOFT_RESET    0xDEAFu

/*
 * PUBLIC API
 ******************************************************************************
 */

/**
 * @brief  Initialise pinmux and I2C peripheral for BMI323.
 *         Call once at startup before any other bmi323_* function.
 * @return 0 on success, -1 on error.
 */
int8_t bmi323_init(void);

/**
 * @brief  Read CHIP_ID register and verify it equals 0x43.
 * @return 0 if CHIP_ID == 0x43, -1 otherwise.
 */
int8_t bmi323_check_id(void);

/**
 * @brief  Issue a soft reset (write 0xDEAF to CMD register) and wait
 *         the required 450 µs POR idle time.
 * @return 0 on success, -1 on I2C error.
 */
int8_t bmi323_reset(void);

/**
 * @brief  Configure accelerometer and gyroscope in low-power duty-cycle mode
 *         at 0.78 Hz.
 * @return 0 on success, -1 on first I2C error.
 */
int8_t bmi323_configure(void);

/**
 * @brief  Read raw accelerometer X, Y, Z counts.
 * @param  x  Output: signed 16-bit raw X count.
 * @param  y  Output: signed 16-bit raw Y count.
 * @param  z  Output: signed 16-bit raw Z count.
 * @return 0 on success, -1 on I2C error.
 */
int8_t bmi323_read_accel(int16_t *x, int16_t *y, int16_t *z);

/**
 * @brief  Read raw gyroscope X, Y, Z counts.
 * @param  x  Output: signed 16-bit raw X count.
 * @param  y  Output: signed 16-bit raw Y count.
 * @param  z  Output: signed 16-bit raw Z count.
 * @return 0 on success, -1 on I2C error.
 */
int8_t bmi323_read_gyro(int16_t *x, int16_t *y, int16_t *z);

/**
 * @brief  Read raw temperature register and print converted value.
 *         T(°C) = (int16_t)raw / 512.0 + 23.0; 0x8000 = invalid.
 * @param  temp  Output: raw 16-bit temperature register value.
 * @return 0 on success, -1 on I2C error.
 */
int8_t bmi323_read_temp(int16_t *temp);

/**
 * @brief  Read one 16-bit register.
 *         Handles the 2-dummy-byte quirk internally.
 * @param  reg  8-bit register address.
 * @param  val  Output: 16-bit register value (little-endian).
 * @return 0 on success, -1 on I2C error.
 */
int8_t bmi323_reg_read(uint8_t reg, uint16_t *val);

/**
 * @brief  Write one 16-bit register.
 * @param  reg  8-bit register address.
 * @param  val  16-bit value to write (little-endian).
 * @return 0 on success, -1 on I2C error.
 */
int8_t bmi323_reg_write(uint8_t reg, uint16_t val);

#endif /* _BMI323_H_ */
