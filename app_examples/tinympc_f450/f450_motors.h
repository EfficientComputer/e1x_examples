#pragma once
#include <stdint.h>

/*
 * F450 motor output: map the four TinyMPC commands to ESC PWM on a PCA9685.
 *
 * The controller emits z[0..3] = per-motor thrust, normalized by hover thrust,
 * delta about hover (Q20): z=0 -> hover, z=+ONE -> +1 hover-thrust. This layer
 * maps each to an ESC pulse width and drives PCA9685 channels 0-3 over I2C.
 *
 * HOVER_US / GAIN_US in f450_motors.c are per-airframe CALIBRATION constants
 * (hover throttle pulse, and us per unit normalized-thrust). The mapping is a
 * deliberately simple linear one; refine it with your motor/prop thrust curve.
 *
 * SAFETY: motors_init arms the ESCs at minimum throttle. Always bench-test with
 * PROPS OFF first. On any abort, call f450_motors_disarm().
 */

/* PCA9685 init + PWM freq + arm all four ESCs at MIN. Returns 0 ok, -1 on I2C error. */
int8_t f450_motors_init(void);

/* Map z[0..3] (Q20 normalized thrust delta) -> pulse us -> PCA9685 ch 0..3. */
void   f450_motors_write(const int32_t z[4]);

/* Failsafe: drive all outputs off (motors stop). */
void   f450_motors_disarm(void);
