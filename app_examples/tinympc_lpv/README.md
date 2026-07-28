# TinyMPC LPV: Fixed-Wing MPC Across a Flight Envelope

This example runs cache-scheduled (LPV) model-predictive control for a fixed-wing
aircraft on E1x. It extends [TinyMPC](https://tinympc.org) from a single
operating point to a whole **flight envelope**, while keeping the online cost
identical.

---

## The idea

Vanilla TinyMPC precomputes one infinite-horizon cache (`Kinf, Pinf, Quu⁻¹,
AmBKt`) from the DARE. This is extremely cheap online, but only valid in a
neighborhood of **one** operating point. That's fine for a quadrotor hovering; a
fixed-wing aircraft, though, flies across a range of airspeeds, climb angles, and
bank angles, and a single linearization goes stale fast.

The fix is classic gain scheduling, applied to TinyMPC's cache: precompute a
**grid of caches**, one per trimmed flight condition, and at runtime **select the
nearest cache** for the live condition, then run the *same* ADMM solve on it.

- Scheduling is a nearest-node array lookup: `~10²` integer ops vs `~10⁵` for a
  solve, so it's free relative to the solve.
- The online solve is unchanged: same fabric kernel, same latency.
- The only cost is the read-only cache LUT (here 21 nodes, ~17 KB in Q20).

This example uses an **Ultra Stick 25e** (a UMN UAV-lab research aircraft) with a
grid of **coordinated-turn** trims over bank × climb (8 states `[u v w, p q r,
phi theta]`, 4 controls `[throttle, elevator, aileron, rudder]`, horizon 15).

---

## What it does

`main.c`:
1. Commands a flight condition (a 30° right coordinated turn).
2. `lpv_select_node` picks the nearest cache node for it.
3. Builds the per-control delta-box about that node's trim control.
4. Solves one ADMM step from a small attitude error, on the fabric kernel and a
   bit-identical scalar reference, and verifies they match (`PASS`/`FAIL`).

The fabric kernel (`tinympc_lpv.c`) is the same Q20 ADMM as the `tinympc`
example, dimensioned for the 8-state aircraft, with a per-control box.

---

## Code structure

- **`tinympc_lpv.c`**: `lpv_select_node` (the schedule) + the `__efficient__`
  fabric ADMM kernel + a scalar reference.
- **`ustick_turn_q20.h`**: the bank × climb cache LUT (21 nodes: `A, B, Kinf,
  Quu⁻¹, AmBKt`), the per-node `(bank, climb)` schedule coordinates, and per-node
  trim controls, all Q20.
- **`main.c`**: schedule, solve, compare.

---

## Build & run

Standard Efficient SDK flow (`fabric` and `scalar` targets); run in sim or flash to the
EVK. Prints the scheduled node, the control delta, and `PASS`/`FAIL`.

---

## Why this example is useful

It shows how to make MPC track a vehicle across its **operating envelope** on
E1x without paying more online: the fabric does the same ~constant-cost ADMM
solve, scheduling just selects which precomputed cache feeds it. This is the
practical recipe for real fixed-wing / variable-condition flight control where a
single LTI model isn't enough. It does so at fixed-point, fabric-accelerated,
low-energy cost. (It covers the quasi-static envelope of trimmed flight
conditions; genuine
non-equilibrium aerobatics are out of scope for any trim-scheduled approach.)
