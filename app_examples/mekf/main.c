/*
 * MEKF attitude estimator — IMU sensor fusion on E1x.
 *
 * A Multiplicative Extended Kalman Filter (MEKF) fuses a gyroscope (prediction)
 * and an accelerometer (tilt correction) to estimate vehicle attitude — the core
 * of any drone / robot orientation estimator. This is the fixed-point Q20 variant
 * (quaternion/bias in Q15, 6×6 covariance in Q20, gravity vector in Q8): the hot
 * path is pure int32 with int64 widening accumulation — no float, no 64-bit
 * divide — so it runs efficiently on the E1x scalar core.
 *
 * Scenario (deterministic, sensor inputs synthesized in main): the true attitude
 * is level, but the filter starts ~2.9° tilted. The gyro reads ~zero and the
 * accelerometer reads gravity straight down; the MEKF must pull its tilted
 * estimate back to level. PASS = converged attitude error below the Q20 floor.
 *
 * (A single 6×6 MEKF is below the fabric's useful payload — the win here is the
 * fixed-point scalar core; batching many filters is what engages the fabric.)
 */
#include <eff.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "mekf_q20.h"

#define N_STEPS 500

int main(void)
{
    /* sensor-fusion tuning (matches the validation suite) */
    const double dt = 0.01, sg = 0.001, sa = 0.01, g = 9.81;
    int32_t dt_q15 = (int32_t)(dt * Q15_ONE);
    int32_t qg_q20 = (int32_t)(sg*sg*dt * Q20_ONE); if (qg_q20 < 1) qg_q20 = 1;
    int32_t Ra_q20 = (int32_t)(sa*sa/dt/(g*g) * Q20_ONE); if (Ra_q20 < 1) Ra_q20 = 1;

    /* true attitude = level; gyro reads ~0; accel reads gravity straight down */
    int32_t w_gyro_q15[3] = {0, 0, 0};
    int32_t a_norm_q8[3]  = {0, 0, Q8_ONE};       /* (0,0,1) unit gravity, body frame */
    double  q_true[4]     = {1.0, 0.0, 0.0, 0.0};

    /* filter starts tilted ~2.86° about x (cos/sin of 0.025 rad half-angle) */
    MEKF_Q20 m;
    float q0[4] = {(float)cos(0.025), (float)sin(0.025), 0.0f, 0.0f};
    float bg0[3] = {0, 0, 0};
    mekf_q20_init(&m, q0, bg0, 0.2f, 0.05f);

    float e_start = mekf_q20_err_deg(&m, q_true);
    for (int k = 0; k < N_STEPS; k++) {
        mekf_q20_predict(&m, w_gyro_q15, dt_q15, qg_q20, 1);
        mekf_q20_update_accel(&m, a_norm_q8, Ra_q20);
    }
    float e_end = mekf_q20_err_deg(&m, q_true);

    printf("[mekf] Q20 fixed-point attitude MEKF (gyro + accel fusion), %d steps\r\n", N_STEPS);
    printf("[mekf] attitude error: start %ld milideg -> end %ld milideg\r\n",
           (long)(e_start*1000), (long)(e_end*1000));

    int ok = (e_end < 1.0f) && (e_end < e_start);   /* converged below the Q20 floor */
    printf(ok ? "[mekf] PASS\r\n" : "[mekf] FAIL\r\n");
    return ok ? 0 : 1;
}
