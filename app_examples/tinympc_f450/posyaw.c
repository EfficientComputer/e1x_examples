/*
 * Fixed-point position/velocity + mag-yaw aiding (see posyaw.h). All int32.
 */
#include "posyaw.h"

#define DT_Q15   655          /* 0.02 s in Q15 */
#define G_Q20    10285974     /* 9.81 m/s² in Q20 */

/* complementary gains (Q15). Baro aids z every step; GPS aids x,y on each fix. */
#define KP_BARO  6554         /* 0.20 */
#define KV_BARO  1311         /* 0.04 */
#define KP_GPS   3277         /* 0.10 */
#define KV_GPS    655         /* 0.02 */

/* ── quaternion (Q15) -> rotation matrix (Q15), body->world ───────────────── */
static void rotmat_q15(const int32_t q[4], int32_t R[9])
{
    int32_t w=q[0], x=q[1], y=q[2], z=q[3];
    int32_t xx=(x*x)>>15, yy=(y*y)>>15, zz=(z*z)>>15;
    int32_t xy=(x*y)>>15, xz=(x*z)>>15, yz=(y*z)>>15;
    int32_t wx=(w*x)>>15, wy=(w*y)>>15, wz=(w*z)>>15;
    R[0]=32768-2*(yy+zz); R[1]=2*(xy-wz);     R[2]=2*(xz+wy);
    R[3]=2*(xy+wz);       R[4]=32768-2*(xx+zz); R[5]=2*(yz-wx);
    R[6]=2*(xz-wy);       R[7]=2*(yz+wx);     R[8]=32768-2*(xx+yy);
}

void posyaw_rotate(const int32_t q[4], const int32_t v[3], int32_t out[3])
{
    int32_t R[9]; rotmat_q15(q, R);
    for (int i = 0; i < 3; i++) {
        int64_t s = (int64_t)R[i*3+0]*v[0] + (int64_t)R[i*3+1]*v[1] + (int64_t)R[i*3+2]*v[2];
        out[i] = (int32_t)(s >> 15);
    }
}

/* ── CORDIC atan2 -> Q15 radians ──────────────────────────────────────────── */
/* atan(2^-i) table in Q15 radians */
static const int32_t ATAN_Q15[15] = {
    25735,15192,8018,4068,2045,1024,512,256,128,64,32,16,8,4,2
};
int32_t posyaw_atan2_q15(int32_t y, int32_t x)
{
    if (x == 0 && y == 0) return 0;
    int32_t ang = 0, xi = x, yi = y;
    /* fold into right half-plane; track the quadrant offset */
    if (xi < 0) {
        int32_t nx=-xi, ny=-yi; xi=nx; yi=ny;
        ang = (y >= 0) ? 102944 : -102944;   /* ±π */
    }
    for (int i = 0; i < 15; i++) {
        int32_t xn, yn;
        if (yi > 0) { xn = xi + (yi>>i); yn = yi - (xi>>i); ang += ATAN_Q15[i]; }
        else        { xn = xi - (yi>>i); yn = yi + (xi>>i); ang -= ATAN_Q15[i]; }
        xi = xn; yi = yn;
    }
    return ang;
}

int32_t posyaw_mag_heading_q15(const int32_t q[4], const int32_t mag_body[3])
{
    int32_t mw[3]; posyaw_rotate(q, mag_body, mw);   /* field in world frame */
    return posyaw_atan2_q15(mw[1], mw[0]);           /* atan2(east, north) */
}

/* ── position/velocity complementary filter ───────────────────────────────── */
void poskf_init(poskf_t *pk)
{
    for (int i = 0; i < 3; i++) { pk->r[i] = 0; pk->v[i] = 0; }
}

void poskf_predict(poskf_t *pk, const int32_t q[4], const int32_t a_body_q20[3])
{
    int32_t aw[3];
    posyaw_rotate(q, a_body_q20, aw);        /* specific force in world (Q20) */
    aw[2] -= G_Q20;                          /* remove gravity -> inertial accel */
    for (int i = 0; i < 3; i++) {
        pk->v[i] += (int32_t)(((int64_t)aw[i]    * DT_Q15) >> 15);
        pk->r[i] += (int32_t)(((int64_t)pk->v[i] * DT_Q15) >> 15);
    }
}

void poskf_update_baro(poskf_t *pk, int32_t baro_z_q20)
{
    int32_t e = baro_z_q20 - pk->r[2];
    pk->r[2] += (int32_t)(((int64_t)KP_BARO * e) >> 15);
    pk->v[2] += (int32_t)(((int64_t)KV_BARO * e) >> 15);
}

void poskf_update_gps(poskf_t *pk, int32_t gx_q20, int32_t gy_q20)
{
    int32_t ex = gx_q20 - pk->r[0], ey = gy_q20 - pk->r[1];
    pk->r[0] += (int32_t)(((int64_t)KP_GPS * ex) >> 15);
    pk->v[0] += (int32_t)(((int64_t)KV_GPS * ex) >> 15);
    pk->r[1] += (int32_t)(((int64_t)KP_GPS * ey) >> 15);
    pk->v[1] += (int32_t)(((int64_t)KV_GPS * ey) >> 15);
}
