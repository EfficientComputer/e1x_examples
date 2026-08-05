/**
 ******************************************************************************
 *
 * @file bmm350.c
 *
 * @brief Bosch BMM350 3-axis magnetometer driver
 *
 * Implements I2C register access for the BMM350 magnetic sensor
 * on the Efficient Computer E1x rev3 EVK.
 *
 * Hardware quirk: BMM350 inserts 2 dummy bytes before payload on every
 * I2C read (datasheet §9.2.3). Callers must request 2 extra bytes and
 * read the payload starting from buf[2].
 *
 ******************************************************************************
 */

/*
 * INCLUDE FILES
 ******************************************************************************
 */
#include "bmm350.h"
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
 * Spin-loop delays.
 * BMM350 POR / reset settling time: ~3 ms.
 * At 50 MHz CPU: 150 000 cycles ≈ 3 ms.
 *
 * Suspend → normal mode startup: ~70 ms.
 * At 50 MHz CPU: 3 500 000 cycles ≈ 70 ms.
 */
#define BMM350_POR_DELAY_CYCLES     150000u
#define BMM350_NORMAL_DELAY_CYCLES  3500000u

/*
 * HELPERS
 ******************************************************************************
 */

/**
 * @brief  Reconstruct a signed 21-bit integer from three raw bytes.
 *         Byte order from the register map: XLSB (bits 7:0), LSB (bits 15:8),
 *         MSB (bits 23:16).  Only bits [20:0] carry data; sign-extend from
 *         bit 20.
 */
static int32_t bmm350_s21(uint8_t xlsb, uint8_t lsb, uint8_t msb)
{
    int32_t raw = ((int32_t)msb << 16) | ((int32_t)lsb << 8) | (int32_t)xlsb;
    /* Mask to 21 bits; bits [23:21] of msb are reserved and may be non-zero */
    raw &= (int32_t)0x001FFFFF;
    /* sign-extend from bit 20 */
    if (raw & (1u << 20))
        raw |= (int32_t)0xFFE00000;
    return raw;
}

/*
 * REGISTER ACCESS
 ******************************************************************************
 */

int8_t bmm350_reg_read(uint8_t reg, uint8_t *val)
{
    /*
     * BMM350 inserts 2 dummy bytes before each payload byte on I2C reads.
     * Request 3 bytes (2 dummy + 1 payload); payload is at buf[2].
     */
    volatile uint8_t buf[4] = {0};
    int8_t rc = eff_i2c_read(I2C_2_1, BMM350_I2C_ADDR, reg, (uint8_t *)buf, 4);
    if (rc != 0) {
        printf("bmm350: read 0x%02x failed\n", (unsigned int)reg);
        return -1;
    }
    *val = buf[2];
    return 0;
}

int8_t bmm350_reg_write(uint8_t reg, uint8_t val)
{
    uint8_t buf[1] = { val };
    int8_t rc = eff_i2c_write(I2C_2_1, BMM350_I2C_ADDR, reg, buf, 1);
    if (rc != 0) {
        printf("bmm350: write 0x%02x failed\n", (unsigned int)reg);
        return -1;
    }
    return 0;
}

/*
 * PUBLIC FUNCTIONS
 ******************************************************************************
 */

int8_t bmm350_init(void)
{
    eff_pinmux_set(PINMUX_2, PINMUX_SPI_I2C1);

    if (eff_i2c_init(I2C_2_1, I2C_SPEED_100K) != 0) {
        printf("bmm350_init: I2C init failed\n");
        return -1;
    }
    eff_gpio_pull_set(GPIO_2, GPIO_PIN_4 | GPIO_PIN_5, EFF_GPIO_PULL_NONE);

    printf("bmm350_init: OK\n");
    return 0;
}

int8_t bmm350_check_id(void)
{
    uint8_t id = 0;
    if (bmm350_reg_read(BMM350_REG_CHIP_ID, &id) != 0)
        return -1;

    if (id != BMM350_CHIP_ID_EXPECTED) {
        printf("bmm350_check_id: unexpected CHIP_ID 0x%02x", (unsigned int)id);
        printf(" (expected 0x%02x)\n", (unsigned int)BMM350_CHIP_ID_EXPECTED);
        return -1;
    }

    printf("bmm350_check_id: CHIP_ID = 0x%02x — BMM350 confirmed\n",
           (unsigned int)id);
    return 0;
}

int8_t bmm350_reset(void)
{
    /* Trigger POR: write 0xB6 to CMD register */
    if (bmm350_reg_write(BMM350_REG_CMD, BMM350_CMD_SOFT_RESET) != 0)
        return -1;

    sleep_ms(4);  /* POR sequence: ~3 ms per datasheet */

    /* Clear CMD register after reset */
    if (bmm350_reg_write(BMM350_REG_CMD, 0x00) != 0)
        return -1;

    sleep_ms(4);

    printf("bmm350_reset: done\n");
    return 0;
}

int8_t bmm350_configure(void)
{
    /*
     * Power off the OTP interface.  After every POR/soft-reset the OTP block
     * stays powered on, which holds cmd_busy=1 indefinitely and prevents the
     * PMU state machine from completing any command.  This single write frees
     * the PMU.  The Bosch SensorAPI issues this implicitly at the end of
     * bmm350_init(); without it no PMU command ever takes effect.
     */
    if (bmm350_reg_write(BMM350_REG_OTP_CMD, BMM350_OTP_CMD_PWR_OFF) != 0)
        return -1;
    sleep_ms(2);

    /*
     * Magnetic reset: BR first then FGR (Bosch SensorAPI order).
     * BR (~12 ms) resets digital filters; FGR (~14 ms) resets the flux gate.
     */
    if (bmm350_reg_write(BMM350_REG_PMU_CMD, BMM350_PMU_CMD_BR) != 0)
        return -1;
    sleep_ms(15);

    if (bmm350_reg_write(BMM350_REG_PMU_CMD, BMM350_PMU_CMD_FGR) != 0)
        return -1;
    sleep_ms(20);

    /* Set ODR = 25 Hz with no averaging, then commit */
    if (bmm350_reg_write(BMM350_REG_PMU_CMD_AGGR_SET,
                         BMM350_AGGR_SET_25HZ) != 0)
        return -1;

    if (bmm350_reg_write(BMM350_REG_PMU_CMD,
                         BMM350_PMU_CMD_UPDATE_OAE) != 0)
        return -1;

    sleep_ms(30);

    /* Enable normal mode */
    if (bmm350_reg_write(BMM350_REG_PMU_CMD, BMM350_PMU_CMD_NORMAL) != 0)
        return -1;

    /* Wait for suspend → normal mode transition (~70 ms per datasheet) */
    sleep_ms(80);

    /* Read PMU status to confirm mode change */
    uint8_t pmu_st = 0;
    bmm350_reg_read(BMM350_REG_PMU_CMD_STATUS_0, &pmu_st);
    printf("bmm350_configure: PMU_STATUS_0 = 0x%02x\n", (unsigned int)pmu_st);

    /* At 25 Hz, first sample takes ~40 ms — wait for data ready */
    sleep_ms(80);

    uint8_t drdy = 0;
    bmm350_reg_read(BMM350_REG_INT_STATUS, &drdy);
    printf("bmm350_configure: INT_STATUS = 0x%02x\n", (unsigned int)drdy);

    printf("bmm350_configure: normal mode, 25 Hz, avg=0\n");
    return 0;
}

int8_t bmm350_read_mag(int32_t *x, int32_t *y, int32_t *z)
{
    /*
     * Burst read starting at MAG_X_XLSB (0x31).
     * Registers auto-increment: 0x31..0x39 = 9 data bytes (XYZ × 3).
     * BMM350 inserts 2 dummy bytes before payload; request 11 total.
     * Layout: buf[0..1]=dummy, buf[2..4]=X, buf[5..7]=Y, buf[8..10]=Z.
     */
    volatile uint8_t buf[12] = {0};
    int8_t rc = eff_i2c_read(I2C_2_1, BMM350_I2C_ADDR,
                              BMM350_REG_MAG_X_XLSB, (uint8_t *)buf, 12);
    if (rc != 0) {
        printf("bmm350: mag burst read failed\n");
        return -1;
    }

    *x = bmm350_s21(buf[2], buf[3], buf[4]);
    *y = bmm350_s21(buf[5], buf[6], buf[7]);
    *z = bmm350_s21(buf[8], buf[9], buf[10]);
    return 0;
}

int8_t bmm350_read_temp(int32_t *temp)
{
    /*
     * Burst read starting at TEMP_XLSB (0x3A).
     * Registers: 0x3A, 0x3B, 0x3C = 3 data bytes.
     * BMM350 inserts 2 dummy bytes before payload; request 5 total.
     * Layout: buf[0..1]=dummy, buf[2..4]=TEMP(XLSB,LSB,MSB).
     */
    volatile uint8_t buf[6] = {0};
    int8_t rc = eff_i2c_read(I2C_2_1, BMM350_I2C_ADDR,
                              BMM350_REG_TEMP_XLSB, (uint8_t *)buf, 6);
    if (rc != 0) {
        printf("bmm350: temp read failed\n");
        return -1;
    }

    *temp = bmm350_s21(buf[2], buf[3], buf[4]);
    return 0;
}
