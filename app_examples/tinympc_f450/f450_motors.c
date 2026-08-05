/*
 * F450 motor output — TinyMPC commands -> ESC PWM via PCA9685 (channels 0-3).
 * See f450_motors.h for the contract and the SAFETY note (props off to bench-test).
 */
#include "f450_motors.h"
#include "pca9685.h"
#include "f450_cache_q20.h"   /* F450Q (Q20 shift) */

#define ADDR      PCA9685_ADDR_DEFAULT
#define FREQ_HZ   400u        /* ESC update rate; 400 Hz keeps actuator latency low.
                               * Use 50u for older analog ESCs that need it.       */
#define MIN_US    1000u       /* ESC idle / arm pulse         */
#define MAX_US    2000u       /* ESC full throttle            */
#define HOVER_US  1400        /* CALIBRATE: pulse at hover throttle for your airframe */
#define GAIN_US   600         /* CALIBRATE: us per unit (Q20) of normalized thrust    */

static uint16_t clamp_us(int v)
{
    if (v < (int)MIN_US) return MIN_US;
    if (v > (int)MAX_US) return MAX_US;
    return (uint16_t)v;
}

int8_t f450_motors_init(void)
{
    if (pca9685_init(ADDR) != 0) return -1;
    if (pca9685_set_freq(ADDR, FREQ_HZ) != 0) return -1;
    for (uint8_t ch = 0; ch < 4; ch++)          /* arm: hold all ESCs at MIN */
        pca9685_set_servo_us(ADDR, ch, MIN_US, FREQ_HZ);
    return 0;
}

void f450_motors_write(const int32_t z[4])
{
    for (uint8_t i = 0; i < 4; i++) {
        /* pulse = HOVER + z_i * GAIN  (z_i is Q20; shift down by F450Q) */
        int us = HOVER_US + (int)(((int64_t)z[i] * GAIN_US) >> F450Q);
        pca9685_set_servo_us(ADDR, i, clamp_us(us), FREQ_HZ);
    }
}

void f450_motors_disarm(void)
{
    pca9685_all_off(ADDR);
}
