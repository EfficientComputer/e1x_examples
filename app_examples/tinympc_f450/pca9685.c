/*
 * pca9685.c — NXP PCA9685 16-channel 12-bit PWM controller driver
 *
 * All 16 channels share a single PWM frequency set by the PRE_SCALE register.
 * Auto-increment is always enabled so multi-byte channel writes are one
 * I2C transaction (4 bytes: ON_L, ON_H, OFF_L, OFF_H).
 *
 * Oscillator: internal 25 MHz (no external crystal needed).
 */

#include "pca9685.h"
#include <eff.h>
#include <eff/drivers/i2c.h>
#include <stdio.h>
#include <stdint.h>

/* PCA9685 I2C bus — DIGIO044=SCL (GPIO_4, PIN_4), DIGIO045=SDA (GPIO_4, PIN_5) */
extern eff_i2c_t *I2C_4_1;

/* ── Register addresses ───────────────────────────────────────────────────── */
#define REG_MODE1           0x00u
#define REG_MODE2           0x01u
#define REG_LED0_ON_L       0x06u   /* base of per-channel regs; ch n = +4*n */
#define REG_ALL_LED_ON_L    0xFAu   /* applies to all 16 channels at once    */
#define REG_PRE_SCALE       0xFEu

/* MODE1 bits */
#define MODE1_RESTART       0x80u
#define MODE1_AI            0x20u   /* auto-increment                        */
#define MODE1_SLEEP         0x10u

/* MODE2 bits */
#define MODE2_OUTDRV        0x04u   /* 1 = totem-pole (needed for servo load) */

/* OFF_H full-off bit — overrides PWM, drives output permanently low */
#define FULL_OFF_BIT        0x10u

/* ~600 µs busy-wait at ~50 MHz e1x clock (datasheet: wait > 500 µs after
 * clearing SLEEP before sending RESTART) */
#define WAKE_CYCLES         30000u

/* ── Helpers ──────────────────────────────────────────────────────────────── */
static int8_t reg_write(uint8_t addr, uint8_t reg, uint8_t val)
{
    /* Use eff_i2c_write_raw directly instead of eff_i2c_write.
     * eff_i2c_write calls eff_i2c_tx_finish which issues ATCIIC100_IIC_RST.
     * That reset fires a spurious INTCMPL and sets _isr_flag=1; the next
     * i2c_wait_int then returns immediately on the stale flag, racing the
     * in-flight transaction and reading a wrong STATUS.ACK.
     * eff_i2c_write_raw sends the full transaction (START+addr+reg+data+STOP)
     * and leaves the bus idle — no IIC_RST, no stale flag. */
    volatile uint8_t data = val;
    return eff_i2c_write_raw(I2C_4_1, addr, reg, I2C_REG, &data, 1);
}

/*
 * prescale = round(osc / (4096 * freq)) - 1
 *
 * Integer rounding via (num + denom/2) / denom:
 *   = (osc + 2048*freq) / (4096*freq) - 1
 * Clamped to [3, 255] per datasheet limits (3 → ~1526 Hz, 255 → ~24 Hz).
 */
static uint8_t calc_prescale(uint16_t freq_hz)
{
    uint32_t val = (PCA9685_OSC_HZ + 2048UL * freq_hz) / (4096UL * freq_hz);
    if (val < 4)   val = 4;   /* clamp: prescale = val-1 >= 3 */
    if (val > 256) val = 256;
    return (uint8_t)(val - 1u);
}

/* ── Public API ───────────────────────────────────────────────────────────── */

int8_t pca9685_init(uint8_t addr)
{
    eff_pinmux_set(PINMUX_4, PINMUX_SPI_I2C1);
    if (eff_i2c_init(I2C_4_1, I2C_SPEED_100K) != 0) {
        printf("pca9685_init: I2C init failed\n");
        return -1;
    }
    eff_gpio_pull_set(GPIO_4, GPIO_PIN_4 | GPIO_PIN_5, EFF_GPIO_PULL_NONE);

    if (reg_write(addr, REG_MODE1, MODE1_SLEEP) != 0) {
        printf("pca9685_init: no ACK at 0x%02X\n", (unsigned int)addr);
        return -1;
    }

    /* Totem-pole output drivers (required to source/sink servo current) */
    if (reg_write(addr, REG_MODE2, MODE2_OUTDRV) != 0) return -1;

    /* Default to 50 Hz (servo standard) */
    if (pca9685_set_freq(addr, 50) != 0) return -1;

    /* All channels off until explicitly set */
    if (pca9685_all_off(addr) != 0) return -1;

    printf("pca9685_init: OK (addr=0x%02X)\n", (unsigned int)addr);
    return 0;
}

int8_t pca9685_set_freq(uint8_t addr, uint16_t freq_hz)
{
    uint8_t prescale = calc_prescale(freq_hz);

    /* PRE_SCALE can only be changed while the oscillator is sleeping */
    if (reg_write(addr, REG_MODE1, MODE1_SLEEP)            != 0) return -1;
    if (reg_write(addr, REG_PRE_SCALE, prescale)           != 0) return -1;

    /* Wake: AI on, sleep off */
    if (reg_write(addr, REG_MODE1, MODE1_AI)               != 0) return -1;

    /* Wait > 500 µs for oscillator to stabilise */
    for (volatile uint32_t i = 0; i < WAKE_CYCLES; i++) {}

    /* Restart all channels so they begin cleanly from tick 0 */
    if (reg_write(addr, REG_MODE1, MODE1_RESTART | MODE1_AI) != 0) return -1;

    return 0;
}

int8_t pca9685_set_channel(uint8_t addr, uint8_t ch, uint16_t on, uint16_t off)
{
    if (ch > 15u || on > 4095u || off > 4095u) return -1;

    uint8_t buf[4] = {
        (uint8_t)(on  & 0xFFu),
        (uint8_t)(on  >> 8u),
        (uint8_t)(off & 0xFFu),
        (uint8_t)(off >> 8u),
    };

    /* Auto-increment: one 4-byte write covers ON_L, ON_H, OFF_L, OFF_H */
    return eff_i2c_write(I2C_4_1, addr, REG_LED0_ON_L + (uint8_t)(ch * 4u),
                         buf, 4);
}

int8_t pca9685_set_servo_us(uint8_t addr, uint8_t ch,
                             uint16_t us, uint16_t freq_hz)
{
    /*
     * off_counts = us * freq_hz * 4096 / 1_000_000
     * Rounding: add 500_000 before dividing.
     * Max safe product (servo use): 2500 µs × 400 Hz × 4096 = ~4.1e9 — fits
     * in uint32_t only just; clamp to 4095 if the result overflows 12 bits.
     */
    uint32_t off = ((uint32_t)us * (uint32_t)freq_hz * 4096u + 500000u)
                   / 1000000u;
    if (off > 4095u) off = 4095u;
    return pca9685_set_channel(addr, ch, 0u, (uint16_t)off);
}

int8_t pca9685_all_off(uint8_t addr)
{
    /*
     * Write to ALL_LED registers (0xFA–0xFD) in one transaction.
     * Setting bit 4 of OFF_H (FULL_OFF) overrides PWM on every channel.
     */
    uint8_t buf[4] = {0x00u, 0x00u, 0x00u, FULL_OFF_BIT};
    return eff_i2c_write(I2C_4_1, addr, REG_ALL_LED_ON_L, buf, 4);
}

