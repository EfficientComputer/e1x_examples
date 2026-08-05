/**
 ******************************************************************************
 *
 * @file bmm350.h
 *
 * @brief Bosch BMM350 3-axis magnetometer driver — public API
 *
 * BMM350 geomagnetic sensor
 * Interface: I2C (8-bit register addresses, 8-bit register values)
 * I2C address: 0x14 (ADSEL=GND) or 0x15 (ADSEL=VDDIO)
 * Host: Efficient Computer E1x rev3 EVK
 *
 * NOTE: BMM350 inserts 2 dummy bytes before payload on every I2C read
 * (datasheet §9.2.3).  Reading any register requires n+2 bytes total:
 * 2 dummy + n data.
 *
 * NOTE: Data output is raw uncompensated 21-bit signed counts.  Full
 * offset/gain/temperature compensation requires the Bosch API OTP
 * coefficients and is not implemented here.
 *
 ******************************************************************************
 */

#ifndef _BMM350_H_
#define _BMM350_H_

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

/* 7-bit I2C address with ADSEL pin tied to GND (bit0 = 0).
 * Change to 0x15u if ADSEL is tied to VDDIO. */
#define BMM350_I2C_ADDR          0x14u

/*
 * MACROS — register addresses
 ******************************************************************************
 */
#define BMM350_REG_CHIP_ID           0x00u
#define BMM350_REG_ERR_REG           0x02u
#define BMM350_REG_PMU_CMD_AGGR_SET  0x04u  /* ODR and averaging config  */
#define BMM350_REG_PMU_CMD_AXIS_EN   0x05u  /* axis enable               */
#define BMM350_REG_PMU_CMD           0x06u  /* power mode command        */
#define BMM350_REG_PMU_CMD_STATUS_0  0x07u  /* PMU status (busy flag)    */
#define BMM350_REG_INT_CTRL          0x2Eu  /* interrupt control; bit7=drdy_data_reg_en */
#define BMM350_REG_INT_STATUS        0x30u  /* data-ready interrupt status; bit2=drdy_data */
#define BMM350_REG_MAG_X_XLSB        0x31u  /* mag X extreme LSB (burst start) */
#define BMM350_REG_TEMP_XLSB         0x3Au  /* temperature extreme LSB   */
#define BMM350_REG_OTP_CMD           0x50u  /* OTP command register      */
#define BMM350_REG_CMD               0x7Eu  /* command register          */

/*
 * MACROS — configuration values
 ******************************************************************************
 */

#define BMM350_CHIP_ID_EXPECTED   0x33u

#define BMM350_CMD_SOFT_RESET     0xB6u   /* triggers POR; follow with 0x00 */
#define BMM350_OTP_CMD_PWR_OFF    0x80u   /* power off OTP interface — must be written
                                            after soft-reset before any PMU command;
                                            while OTP is on cmd_busy never clears */

/* PMU_CMD.pmu_cmd values */
#define BMM350_PMU_CMD_SUSPEND    0x00u
#define BMM350_PMU_CMD_NORMAL     0x01u
#define BMM350_PMU_CMD_UPDATE_OAE 0x02u   /* update ODR/avg after writing AGGR_SET */
#define BMM350_PMU_CMD_FGR        0x05u   /* Flux Gate Reset — required after POR before NORMAL */
#define BMM350_PMU_CMD_BR         0x07u   /* Bit Reset — recommended after POR */

/* PMU_CMD_AXIS_EN: enable axes (en_x=bit0, en_y=bit1, en_z=bit2).
 * Reset value is 0x07 (all enabled), so an explicit write is not required
 * after power-on or soft reset. */
#define BMM350_PMU_CMD_AXIS_EN_ALL 0x07u

/* PMU_CMD_AGGR_SET ODR codes (bits[3:0]):
 * 0x01=400Hz, 0x02=200Hz, 0x03=100Hz, 0x04=50Hz,
 * 0x05=25Hz, 0x06=12.5Hz, 0x07=6.25Hz, 0x08=3.125Hz, 0x09=1.5625Hz
 * bits[7:4] = avg: 0x00 = no averaging */
#define BMM350_AGGR_SET_25HZ  0x05u

/*
 * PUBLIC API
 ******************************************************************************
 */

/**
 * @brief  Initialise pinmux and I2C peripheral for BMM350.
 *         Call once at startup before any other bmm350_* function.
 * @return 0 on success, -1 on error.
 */
int8_t bmm350_init(void);

/**
 * @brief  Read CHIP_ID register and verify it equals 0x33.
 * @return 0 if CHIP_ID == 0x33, -1 otherwise.
 */
int8_t bmm350_check_id(void);

/**
 * @brief  Issue a soft reset (write 0xB6 then 0x00 to CMD register)
 *         and wait for the POR sequence to complete (~3 ms each).
 * @return 0 on success, -1 on I2C error.
 */
int8_t bmm350_reset(void);

/**
 * @brief  Configure BMM350 in normal mode at 25 Hz.
 *         Sets PMU_CMD_AGGR_SET = 0x05 (25 Hz, no averaging),
 *         issues ODR/avg update, then enables normal mode.
 *         Blocks for suspend→normal transition + first data ready.
 * @return 0 on success, -1 on I2C error.
 */
int8_t bmm350_configure(void);

/**
 * @brief  Burst-read all three magnetic-field axes as raw 21-bit signed counts.
 *         Reads 9 data bytes (+ 2 dummy) from registers 0x31–0x39 in one
 *         transaction.  Output is uncompensated raw data.
 * @param  x  Output: signed 21-bit raw X count (sign-extended to int32_t).
 * @param  y  Output: signed 21-bit raw Y count.
 * @param  z  Output: signed 21-bit raw Z count.
 * @return 0 on success, -1 on I2C error.
 */
int8_t bmm350_read_mag(int32_t *x, int32_t *y, int32_t *z);

/**
 * @brief  Read temperature register as raw 21-bit signed count.
 *         Reads 3 data bytes (+ 2 dummy) from registers 0x3A–0x3C.
 * @param  temp  Output: signed 21-bit raw temperature count.
 * @return 0 on success, -1 on I2C error.
 */
int8_t bmm350_read_temp(int32_t *temp);

/**
 * @brief  Read one 8-bit register.
 *         Handles the 2-dummy-byte quirk internally.
 * @param  reg  8-bit register address.
 * @param  val  Output: 8-bit register value.
 * @return 0 on success, -1 on I2C error.
 */
int8_t bmm350_reg_read(uint8_t reg, uint8_t *val);

/**
 * @brief  Write one 8-bit register.
 * @param  reg  8-bit register address.
 * @param  val  8-bit value to write.
 * @return 0 on success, -1 on I2C error.
 */
int8_t bmm350_reg_write(uint8_t reg, uint8_t val);

#endif /* _BMM350_H_ */
