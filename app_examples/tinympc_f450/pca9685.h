#pragma once
#include <stdint.h>

/* I2C address — set by A0-A5 pins; all-low (default) gives 0x40 */
#define PCA9685_ADDR_DEFAULT    0x40u

/* Internal oscillator frequency (Hz).  Datasheet spec: 25 MHz ± 1 % */
#define PCA9685_OSC_HZ          25000000UL

/* Full 12-bit range */
#define PCA9685_COUNTS          4096u

/* Typical servo pulse widths (µs) */
#define PCA9685_SERVO_MIN_US    1000u   /* full deflection one way  */
#define PCA9685_SERVO_MID_US    1500u   /* neutral / centre         */
#define PCA9685_SERVO_MAX_US    2000u   /* full deflection other way */

/*
 * pca9685_init — wake oscillator, enable auto-increment, set 50 Hz default.
 * addr: I2C address (use PCA9685_ADDR_DEFAULT for the standard breakout).
 * Returns 0 on success, -1 on I2C error.
 */
int8_t pca9685_init(uint8_t addr);

/*
 * pca9685_set_freq — set PWM output frequency.
 * Valid range: 24–1526 Hz.  The chip is briefly slept to update PRE_SCALE.
 * Changing frequency while channels are active will glitch outputs momentarily.
 */
int8_t pca9685_set_freq(uint8_t addr, uint16_t freq_hz);

/*
 * pca9685_set_channel — set raw 12-bit ON/OFF tick counts for one channel.
 * ch: 0–15.  on, off: 0–4095 (tick within the 4096-tick period).
 * The output goes high at tick `on` and low at tick `off`.
 * For servo use: set on = 0, off = counts from pca9685_us_to_counts().
 */
int8_t pca9685_set_channel(uint8_t addr, uint8_t ch, uint16_t on, uint16_t off);

/*
 * pca9685_set_servo_us — set servo pulse width in microseconds.
 * freq_hz must match the value last set with pca9685_set_freq().
 * Intended for servo use at 50–400 Hz; higher frequencies with large us values
 * will be clamped to 4095 counts.
 */
int8_t pca9685_set_servo_us(uint8_t addr, uint8_t ch,
                             uint16_t us, uint16_t freq_hz);

/*
 * pca9685_all_off — drive all 16 channels to the full-OFF state immediately.
 * Uses the ALL_LED registers so this is a single 4-byte I2C write.
 */
int8_t pca9685_all_off(uint8_t addr);
