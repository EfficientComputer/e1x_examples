#ifndef MEKF_Q20_H
#define MEKF_Q20_H
#include <stdint.h>

/* Q20 mixed-format MEKF.
 * q, b_g:  Q15 (×32768)
 * P:       Q20 (×1048576) — 32× more precision than Q15
 * h, a:    Q8  (×256)     — unit gravity vector; keeps H×P in int32
 *
 * Key improvement over Q15: cross-covariance P[att,bias] builds at
 * 52 counts/step (vs 1-2 in Q15), so bias Kalman gain is non-zero
 * from step 1 rather than step 170. */

#define Q15_ONE  32768
#define Q20_ONE  1048576
#define Q8_ONE   256

typedef struct {
    int32_t q[4];   /* quaternion, Q15 */
    int32_t b_g[3]; /* gyro bias, Q15 (1.0 = 1 rad/s) */
    int32_t P[36];  /* 6×6 covariance, Q20 */
} MEKF_Q20;

void mekf_q20_init(MEKF_Q20 *m,
                   const float q0[4], const float bg0[3],
                   float sigma_q0, float sigma_bg0);

/* All args pre-converted to their fixed formats:
 *   w_gyro_q15: rad/s in Q15
 *   dt_q15:     seconds in Q15
 *   qg_q20:     sigma_g^2 * dt in Q20, min 1
 *   qbg_q20:    sigma_bg^2 * dt in Q20, min 1 */
void mekf_q20_predict(MEKF_Q20 *m,
                      const int32_t w_gyro_q15[3],
                      int32_t dt_q15,
                      int32_t qg_q20,
                      int32_t qbg_q20);

/* a_norm_q8: a_meas/g in Q8 (unit gravity vector in body frame)
 * Ra_q20:   sigma_a^2/(dt*g^2) in Q20, min 1 */
void mekf_q20_update_accel(MEKF_Q20 *m,
                            const int32_t a_norm_q8[3],
                            int32_t Ra_q20);

float mekf_q20_err_deg(const MEKF_Q20 *m, const double q_ref[4]);
void  mekf_q20_get_bg(const MEKF_Q20 *m, float bg[3]);

#endif
