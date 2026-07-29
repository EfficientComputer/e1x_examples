# QR Decomposition (QR)

This example computes a **QR decomposition** of a matrix on the Electron E1x general-purpose processor, using the Gram-Schmidt process in fixed-point arithmetic. It factors a matrix `A` into an orthogonal matrix `Q` and an upper-triangular matrix `R` such that `A = Q * R`. It is a compact linear-algebra kernel used to show how a column-by-column, dependency-carrying computation maps onto the Fabric architecture.

---

## Overview

### What is a QR Decomposition?

A QR decomposition factors a matrix `A` into an orthogonal matrix `Q` and an upper-triangular matrix `R`, so that `A = Q * R`. The columns of `Q` form an orthonormal basis, and `R` records how the original columns of `A` are expressed in that basis. This example computes it with the Gram-Schmidt process, which builds `Q` one column at a time: each new column of `A` has its projections onto the earlier `Q` columns subtracted off, and the remainder is normalized.

This example works in fixed-point arithmetic: real values are represented as integers scaled by a fixed number of fractional bits, and each multiply is rescaled to keep the result in the same format. The normalization step uses an integer square root. Fixed-point avoids the cost of floating-point hardware, which suits energy-constrained devices.

### Mathematical Definition

For each column `k` of `A`, the Gram-Schmidt process computes:

    R[k][k] = || v_k ||                       (the L2 norm of the working column)

    Q[:,k]  = v_k / R[k][k]                   (normalize to a unit column)

    R[k][j] = Q[:,k] . A[:,j]                 for each later column j

    v_j     = v_j - R[k][j] * Q[:,k]          (subtract the projection from column j)

Where:

- `A` is the input matrix (`n x n` here)
- `Q` is the orthogonal output whose columns are orthonormal
- `R` is the upper-triangular output, satisfying `A = Q * R`

Each column is processed in turn, and computing it depends on the columns finished before it.

---

## Why This Kernel Matters

QR decomposition is one of the most widely used operations in numerical computing:

- **Least-squares fitting**, where QR solves overdetermined systems in a numerically stable way
- **Linear solvers and eigenvalue methods**, where repeated QR steps drive iterative algorithms
- **Signal processing**, where orthogonalization appears in beamforming and adaptive filtering
- **Robotics and graphics**, where orthonormal bases and orientations are extracted from measured data

Because it mixes dot-product reductions, a square root, a division, and projection updates with strong dependencies between columns, it is a good measure of how well an architecture pipelines a structured, dependency-carrying computation.

---

## Why Efficient Hardware Performs Well

The Electron E1x runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

This fits fixed-point Gram-Schmidt QR decomposition well:

- The multiply-accumulate reductions that form the norms and dot products run **in parallel** across the Fabric
- Once a `Q` column is fixed, the projection and update of the later columns are independent and are processed **several at a time**
- Norms, dot products, and the working column stay **resident on the Fabric** instead of spilling to memory
- Fixed-point arithmetic keeps the datapath **simple**, so more of the Fabric is spent on useful work

The result is high arithmetic throughput at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## Configurable Parameters

These definitions in the app headers control the benchmark. Change them to resize the problem or re-run it.

| Definition         | Default | Effect                                                                                                                                                                                                      |
| ------------------ | ------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `NUM_ITERATIONS`   | `1`     | This is how many times the kernel runs, set in `main.c`. Increase it to average out measurement noise when benchmarking.                                                                                            |
| `INPUT_SIZE`       | `16`    | This is the matrix dimension `n` (the input is `n x n`), defined in `qr_decomp.h`. This sets the problem size: a larger value means a larger matrix and more work per run. The input is generated in `qr_decomp.c`. |
| `FIXED_POINT_BITS` | `10`    | This is the number of fractional bits in the fixed-point format, defined in `qr.h`. It sets the trade-off between numeric range and precision. Changing it changes the scaling used throughout the kernel.          |
