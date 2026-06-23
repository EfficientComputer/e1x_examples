# MEKF — IMU Attitude Estimator (Sensor Fusion)

A **Multiplicative Extended Kalman Filter (MEKF)** for vehicle attitude
estimation on E1x — the orientation estimator at the heart of any drone, robot,
or AHRS. It fuses a **gyroscope** (prediction) with an **accelerometer** (gravity
/ tilt correction) to track a quaternion attitude and gyro bias.

> This is sensor *fusion* (an algorithm that consumes IMU data), not a sensor
> chip driver — the IMU inputs here are synthesized in `main.c`. For chip drivers
> see `../../sensor_examples/`. In a real system you'd feed it readings from an
> IMU (e.g. via one of those drivers).

---

## 1. Overview

The MEKF estimates orientation as a nominal quaternion plus a small error state
`[δθ, δb_g]` (attitude error + gyro bias) with a 6×6 covariance:

- **predict** — integrate the gyro to advance the quaternion; propagate covariance.
- **update (accel)** — the accelerometer measures gravity's direction in the body
  frame; the filter corrects tilt (roll/pitch) and learns the gyro bias.

This is the fixed-point **Q20** variant tuned for the E1x:
- quaternion & bias in **Q15**, the 6×6 covariance in **Q20**, gravity vector in
  **Q8** — a mixed format that keeps the cross-covariance terms precise enough to
  estimate bias (a single global Q15 scale is too coarse for the bias terms).
- the hot path (`predict`, `update_accel`) is **pure int32** with int64 widening
  accumulation — **no float, no 64-bit divide**.

---

## 2. What it does

`main.c` runs a deterministic convergence test: the true attitude is **level**,
but the filter is initialized **~2.9° tilted**. The synthesized gyro reads ~zero
and the accelerometer reads gravity straight down. Over 500 steps the MEKF must
pull its estimate back to level.

`PASS` = the attitude error converges below the Q20 precision floor (< 1°).

---

## 3. Code structure

- **`mekf_q20.c` / `mekf_q20.h`** — the fixed-point MEKF: `mekf_q20_init`,
  `mekf_q20_predict`, `mekf_q20_update_accel`, plus float read-outs
  (`mekf_q20_err_deg`, `mekf_q20_get_bg`) used only for reporting.
- **`main.c`** — synthesizes the IMU inputs, runs predict/update, checks
  convergence (`PASS`/`FAIL`).

---

## 4. Build & run

Standard EFF SDK flow (`scalar` target); run in sim or flash to the EVK. Prints
the start/end attitude error and `PASS`/`FAIL`.

---

## 5. Notes on the E1x

A single 6×6 MEKF update is small — below the fabric's useful payload, so the win
here is the **efficient fixed-point scalar core** (no soft-float traps), not
fabric parallelism. To engage the fabric you batch *many* independent estimators
as block-diagonal matrix multiplies (e.g. a swarm, or particle/sigma-point sets),
where the wide work fills the array. This example is the single-filter building
block.
