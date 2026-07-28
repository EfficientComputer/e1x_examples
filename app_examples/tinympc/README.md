# TinyMPC: Quadrotor Model-Predictive Control

This example is a real-time **model-predictive controller (MPC)** for a
Crazyflie-class quadrotor, running on E1x. This is
[TinyMPC](https://tinympc.org), reimplemented in bare-metal C and fixed-point for
the fabric. It is ADMM-based MPC with a precomputed infinite-horizon cache.

---

## Overview

MPC computes control inputs by solving a constrained optimization problem at
every timestep: minimize a quadratic cost over a short horizon subject to the
system dynamics and actuator limits, then apply the first control. It's the
state-of-the-art for agile robots (drones, legged robots, manipulators) because
it respects constraints and looks ahead, but the per-step solve has historically
been too expensive for small, power-constrained flight controllers.

TinyMPC makes it cheap by moving all the heavy linear algebra **offline**: from
the system `(A, B, Q, R, ρ)` it precomputes the infinite-horizon LQR cache
(`Kinf, Pinf, Quu⁻¹, AmBKt`) via the discrete Riccati equation. **Online**, each
step is just a fixed sequence of matrix-vector sweeps (ADMM) with **no matrix
factorization**, which is exactly what maps well onto the E1x fabric.

This example solves the canonical TinyMPC quadrotor problem:
- **12 states** `[x y z, roll pitch yaw, vx vy vz, wx wy wz]`
- **4 controls** (motor thrust deltas), box-constrained
- **horizon N = 10**, **10 ADMM iterations**

---

## Why this maps well to E1x

- The whole multi-iteration solve is one `__efficient__` fabric dispatch.
- All arithmetic is **int32 Q20** fixed point (`int32 × int32 → int64` widening
  accumulate, then `>> 20`), with **no floating point, no 64-bit divide**.
- The cache is read-only data; the online solve is constant-matrix matvecs.

Each ADMM iteration runs:
```
backward:   r = -ρ(z - y);  d = Quu⁻¹(Bᵀp + r);  p = AmBKtᵀp - Kinfᵀr
forward:    u = -Kinf·x - d;  x⁺ = A·x + B·u
slack+dual: z = clip(u + y, [umin, umax]);  y += u - z
```

---

## Code structure

- **`tinympc.c`**: the `__efficient__` fabric ADMM kernel
  (`tinympc_solve_fabric`) plus a bit-identical scalar reference
  (`tinympc_solve_scalar`) for the correctness check.
- **`quad_cache_q20.h`**: the precomputed Crazyflie cache (`A, B, Kinf, Quu⁻¹,
  AmBKt`) and box limits, in Q20.
- **`main.c`**: sets a hover error state, runs both solvers, prints the first
  control and verifies the fabric result matches the reference (`PASS`/`FAIL`).

---

## Build & run

Build with the standard EFF SDK flow (`fabric` and `scalar` targets). Run in the
simulator or flash to the EVK; the app prints the computed control and a
`PASS`/`FAIL` over the UART console.

---

## Correctness

The fabric kernel and the scalar reference perform the identical Q20 arithmetic;
`__effcc_parallel` only parallelizes independent matvec outputs, so the integer
results are **bit-identical**. `PASS` means the fabric solve matches the
reference exactly. (The underlying fixed-point solver is validated against the
floating-point TinyMPC reference elsewhere to sub-millimeter tracking on the
official Crazyflie benchmark.)

---

## Why this example is useful

It demonstrates a **complete, constrained, real-time optimal controller** running
in fixed point on the fabric, not a microbenchmark. The same kernel sustains
real-time rates for quadrotor flight control at a fraction of the energy of an
FPU-based microcontroller, which is exactly the regime where E1x's efficiency
matters: battery-powered, always-on robotics.
