/**
 ******************************************************************************
 *
 * @file ens220.h
 *
 * @brief ScioSense ENS220 Barometric Pressure and Temperature Sensor — public API
 *
 * ENS220 barometric pressure + temperature sensor
 * Interface: I2C with 8-bit register addressing (eff_i2c_read / eff_i2c_write)
 * I2C address: 0x20 (fixed)
 * Host: Efficient Computer E1x rev3 EVK
 *
 ******************************************************************************
 */

#ifndef _ENS220_H_
#define _ENS220_H_

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
#define ENS220_I2C_ADDR      0x20u

/*
 * MACROS — registers
 ******************************************************************************
 */
#define ENS220_REG_PART_ID   0x00u   /* 16-bit LE; expected 0x0321 */
#define ENS220_REG_MODE_CFG  0x06u   /* HP=bit7, START=bit4, MEAS_T=bit1, MEAS_P=bit0, RESET=bit3 */
#define ENS220_REG_MEAS_CFG  0x07u   /* P conversion time / PT rate */
#define ENS220_REG_STBY_CFG  0x08u   /* 0x01 = one-shot */
#define ENS220_REG_OVS_CFG   0x09u   /* oversampling; 0x00 = 1x */
#define ENS220_REG_DATA_STAT 0x14u   /* bit0=TR (T ready), bit1=PR (P ready) */
#define ENS220_REG_PRESS_OUT 0x17u   /* 24-bit unsigned LE; units = 1/64 Pa */
#define ENS220_REG_TEMP_OUT  0x1Au   /* 16-bit unsigned LE; units = 1/128 K */

/*
 * MACROS — MODE_CFG values
 ******************************************************************************
 */
#define ENS220_MODE_RESET    0x08u   /* soft reset; restores idle/reset state */
#define ENS220_MODE_MEASURE  0x93u   /* HP=1, START=1, MEAS_T=1, MEAS_P=1 */

/*
 * MACROS — STBY_CFG
 ******************************************************************************
 */
#define ENS220_STBY_ONESHOT  0x01u   /* auto-idle after 1 measurement */

/*
 * MACROS — expected PART_ID
 ******************************************************************************
 */
#define ENS220_PART_ID_EXPECTED  0x0321u

/*
 * MACROS — spin-loop delay counts (50 MHz CPU clock)
 ******************************************************************************
 */
#define ENS220_RESET_CYCLES   50000u   /* ~1 ms */
#define ENS220_MEAS_CYCLES   500000u   /* ~10 ms (P+T conversion margin) */
#define ENS220_POLL_CYCLES    10000u   /* ~200 µs between DATA_STAT polls */

/*
 * PUBLIC API
 ******************************************************************************
 */

/**
 * @brief  Initialise pinmux and I2C peripheral for ENS220.
 *         Performs a soft reset and configures one-shot standby mode.
 *         Call once at startup before any other ens220_* function.
 * @return 0 on success, -1 on error.
 */
int8_t ens220_init(void);

/**
 * @brief  Read and verify the ENS220 PART_ID.
 *         Expected: 0x0321.
 * @return 0 if ID matches, -1 otherwise.
 */
int8_t ens220_read_id(void);

/**
 * @brief  Issue a soft reset and wait ≥1 ms for the sensor to return to idle.
 * @return 0 on success, -1 on error.
 */
int8_t ens220_reset(void);

/**
 * @brief  Trigger one one-shot measurement and return raw counts.
 *         Polls DATA_STAT until both T and P are ready (up to ~10 ms).
 * @param[out] raw_p  24-bit pressure raw count (units: 1/64 Pa)
 * @param[out] raw_t  16-bit temperature raw count (units: 1/128 K)
 * @return 0 on success, -1 on timeout or error.
 */
int8_t ens220_read_raw(uint32_t *raw_p, uint16_t *raw_t);

/**
 * @brief  Trigger one measurement, print result as integer hPa and °C.
 *           ens220: raw_p=<n> raw_t=<n> P=<n>.<nn> hPa T=<n>.<n> C
 * @return 0 on success, -1 on timeout or error.
 */
int8_t ens220_measure(void);

#endif /* _ENS220_H_ */
