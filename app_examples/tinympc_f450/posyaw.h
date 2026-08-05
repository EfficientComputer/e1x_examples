#pragma once
#include <stdint.h>

/*
 * Fixed-point position/velocity + mag-yaw aiding for the F450 stack — the staged
 * follow-up to the inner loop. All int32 (no float / no soft-float):
 *
 *  - poskf: an inertial + GPS/baro COMPLEMENTARY filter for world-frame position
 *    and velocity. Complementary (integer gains, no covariance matrix) rather than
 *    an EKF, so there is nothing to mis-scale — robust in fixed point.
 *  - mag heading via a CORDIC atan2 (integer), blended into yaw complementarily.
 *
 * Conventions: world z-up; quaternion q is body->world in Q15; positions/velocities
 * in Q20 (m, m/s) to match the controller state. dt is fixed at DT_Q15 (50 Hz).
 */

/* Rotate a body-frame vector (Q20) into the world frame by q (Q15). */
void posyaw_rotate(const int32_t q_q15[4], const int32_t v_body_q20[3], int32_t v_world_q20[3]);

/* CORDIC atan2(y, x) -> angle in Q15 radians (range ±π). x,y any common scale. */
int32_t posyaw_atan2_q15(int32_t y, int32_t x);

/* Magnetometer heading (yaw), Q15 rad: rotate body field to world, atan2(E, N). */
int32_t posyaw_mag_heading_q15(const int32_t q_q15[4], const int32_t mag_body[3]);

/* ── Position/velocity complementary filter ─────────────────────────────── */
typedef struct {
    int32_t r[3];   /* world position, Q20 (m)   */
    int32_t v[3];   /* world velocity, Q20 (m/s) */
} poskf_t;

void poskf_init(poskf_t *pk);

/* Predict: integrate specific force (body, Q20 m/s²) rotated to world minus g. */
void poskf_predict(poskf_t *pk, const int32_t q_q15[4], const int32_t a_body_q20[3]);

/* Correct altitude from baro (Q20 m). Runs every step. */
void poskf_update_baro(poskf_t *pk, int32_t baro_z_q20);

/* Correct horizontal position from a GPS fix (Q20 m, world x=north, y=east). */
void poskf_update_gps(poskf_t *pk, int32_t gps_x_q20, int32_t gps_y_q20);
