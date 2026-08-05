# TinyMPC — Hawk's Work F450 Hover Control

A real-time **model-predictive controller** for a **Hawk's Work F450** (DJI
F450-geometry, ~1.2 kg) quadrotor, running on E1x. Same [TinyMPC](https://tinympc.org)
ADMM solver as the Crazyflie example — only the model (mass, inertia, arm length)
differs, so this example ships a **generator** you can retune to your own airframe.

---

## 1. Overview

The controller regulates the quad back to hover subject to per-motor thrust
limits. It uses the standard 12-state hover linearization:

- **12 states** `[x y z, roll pitch yaw, vx vy vz, wx wy wz]`
- **4 controls** — per-motor thrust, normalized by hover thrust, delta about hover
  (`u = 0` at hover), box-constrained to `|u| ≤ 0.9` (≈ zero-to-double thrust)
- **horizon N = 10**, **10 ADMM iterations**, 50 Hz

Near hover the state coupling `A` depends only on gravity and the timestep; the
**airframe (mass / inertia / geometry) enters entirely through `B`**. That is what
`gen_f450_cache.py` builds before solving the DARE.

## 2. The F450 model (`gen_f450_cache.py`)

Physical parameters (documented estimates for a ~1.2 kg F450 on 3S/4S — **edit
these for your airframe** and regenerate):

| param | value | meaning |
|---|---|---|
| `m`  | 1.2 kg | all-up mass |
| `Ixx=Iyy` | 0.012 kg·m² | roll/pitch inertia |
| `Izz` | 0.022 kg·m² | yaw inertia |
| `L` | 0.225 m | centre-to-motor (450 mm diagonal ÷ 2) |
| `cT` | 0.016 m | drag-torque / thrust ratio (yaw authority) |

```bash
python3 gen_f450_cache.py      # -> f450_cache_q20.h   (needs numpy + scipy)
```

It builds the continuous hover dynamics, discretizes at 50 Hz (matrix
exponential), folds the ADMM penalty `ρ` into the DARE, solves for the
infinite-horizon cache `(Kinf, Pinf, Quu⁻¹, AmBKt)`, and quantizes to **Q20**. It
prints the closed-loop spectral radius as a stability check (currently **0.964**,
stable) and the max gain magnitude (fits the Q20 int32 range).

## 3. Why this maps well to E1x

- The whole multi-iteration solve is one `__efficient__` fabric dispatch.
- All arithmetic is **int32 Q20** (`int32 × int32 → int64` widening accumulate,
  `>> 20`) — **no floating point, no 64-bit divide**.
- The cache is read-only; the online solve is constant-matrix matvecs (ADMM):
```
backward:   r = -ρ(z - y);  d = Quu⁻¹(Bᵀp + r);  p = AmBKtᵀp - Kinfᵀr
forward:    u = -Kinf·x - d;  x⁺ = A·x + B·u
slack+dual: z = clip(u + y, [umin, umax]);  y += u - z
```

## 4. Code structure

- **`gen_f450_cache.py`** — the F450 model + DARE + Q20 quantizer → `f450_cache_q20.h`.
- **`f450_cache_q20.h`** — generated cache (`A, B, Kinf, Quu⁻¹, AmBKt`, box), Q20.
- **`tinympc.c` / `.h`** — the generic `__efficient__` ADMM kernel + bit-identical
  scalar reference (shared with the Crazyflie example; dims are the same).
- **`pca9685.c` / `.h`** — I²C driver for the PCA9685 16-channel PWM chip.
- **`f450_motors.c` / `.h`** — maps the four solver outputs to ESC pulse widths and
  drives PCA9685 channels 0–3 (arm / write / failsafe-disarm).
- **`main.c`** — a **wall-safe "in-place" cage demo** in trajectory-*tracking* mode
  (`e = x - xref(t)`): **hover → vertical bounce → yaw spin-in-place → hover**. It
  only commands motions a tight indoor cage + an IMU/mag/baro/GPS suite handle well
  (vertical is baro-precise; yaw doesn't translate), **never commands horizontal
  translation**, and enforces a **hard geofence** (abort if the craft drifts past a
  safe radius — the real-flight guard against GPS wander/disturbance). Each 50 Hz
  step solves on the fabric, applies the first control, propagates the model, and
  checks the geofence; verifies fabric == scalar on the first solve; `PASS` when it
  tracks the routine and returns to hover without a geofence trip.

## 5. Build & run

Standard EFF SDK flow (`fabric` and `scalar` targets). In sim or on the EVK it
prints the routine (phase, altitude vs reference, yaw, horizontal offset) and a
`PASS`/`FAIL` over the UART console. Expected: it tracks the vertical bounce to
~25 mm, spins a full 360° in place, holds ~0 horizontal offset (so the geofence
never trips), and returns to hover.

### Why these maneuvers (indoor / tight-cage note)

On IMU + mag + baro + GPS, **horizontal position is only GPS-grade (~1 m)** — fine
for large outdoor patterns (≥ ~10 m), but a large fraction of a tight (~4×4×3 m)
cage, and unusable for sub-metre indoor trajectories. This demo therefore only
commands what the sensors resolve well and what can't drift into a wall: **vertical
moves** (baro-precise, ~0.1–0.3 m) and **yaw in place** (no translation). The
geofence is the backstop. For precise *horizontal* indoor work (tight figure-8,
small cross, cm hover) add **optical flow + a downward rangefinder** (~$30, its
sweet spot is exactly low-and-slow) or motion capture, and swap the GPS update in
the estimator for the flow velocity update — then the tight maneuvers open up.

## 6. Motor output (PCA9685 → 4 ESCs)

The E1x has no native PWM peripheral, so the four motor commands go out over I²C to
a **PCA9685** 16-channel PWM chip, which generates the ESC signals (channels 0–3).
This is enabled by the `F450_MOTORS` compile definition (set for the E1x targets in
`CMakeLists.txt`); build without it to run the loop as a pure sim with no PWM.

- **Wiring**: PCA9685 on I²C bus `I2C_4_1` (SCL = DIGIO044, SDA = DIGIO045), address
  `0x40`. ESC signal leads on PCA9685 channels 0–3 (motors 1–4).
- **Mapping** (`f450_motors.c`): `pulse_us = HOVER_US + z_i · GAIN_US`, clamped to
  1000–2000 µs, at `FREQ_HZ` (400 Hz default — low latency; drop to 50 Hz for older
  analog ESCs). `HOVER_US` and `GAIN_US` are **per-airframe calibration** (hover
  throttle pulse, and µs per unit of normalized-thrust command); the linear map is a
  starting point — refine it with your motor/prop thrust curve.
- **Arming / failsafe**: `f450_motors_init()` holds all four ESCs at 1000 µs to arm;
  the geofence abort and the end of the routine both call `f450_motors_disarm()`
  (all channels off).

> ⚠️ **SAFETY — bench-test with PROPS OFF.** Flashing the E1x targets *arms the ESCs
> and emits live PWM* that follows the demo trajectory. In this example the PWM is
> driven from the demo's *simulated* state, so it's for verifying the output path
> (scope the four channels, or run ESCs without props). For real flight, feed the
> controller your *estimated* state instead of the internal model, and keep the
> geofence as the failsafe.

## 7. Correctness

`__effcc_parallel` only parallelizes independent matvec outputs, so the fabric and
scalar integer results are **bit-identical** — `PASS` means they matched exactly and
the closed loop returned to hover. The underlying fixed-point solver is validated
against the floating-point TinyMPC reference (sub-mm on the official Crazyflie
benchmark); this example swaps in the F450 cache, whose closed-loop stability is
checked by the generator (spectral radius < 1) and a hover-recovery rollout.
