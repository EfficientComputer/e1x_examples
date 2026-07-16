# Cholesky Decomposition (Cholesky)

This example computes a **Cholesky decomposition** of a symmetric positive-definite matrix on the Electron E1x general-purpose processor, using fixed-point arithmetic. It factors a matrix `A` into a lower-triangular matrix `L` such that `A = L * L^T`. It is a compact linear-algebra kernel used to show how a triangular, dependency-carrying computation maps onto the Fabric architecture.

---

## 1. Overview

### What is a Cholesky Decomposition?

A Cholesky decomposition takes a symmetric positive-definite matrix `A` and factors it into a lower-triangular matrix `L` and its transpose `L^T`, so that `A = L * L^T`. It is the standard, efficient way to factor such a matrix, using about half the work of a general LU decomposition because it exploits the symmetry.

This example works in fixed-point arithmetic: real values are represented as integers scaled by a fixed number of fractional bits, and each multiply is rescaled to keep the result in the same format. Fixed-point avoids the cost of floating-point hardware, which suits energy-constrained devices.

### Mathematical Definition

For a symmetric positive-definite matrix `A`, the decomposition produces a lower-triangular `L` with:

    L[j][j] = sqrt( A[j][j] - Σ L[j][k]^2 )              for k = 0 .. j-1

    L[i][j] = ( A[i][j] - Σ L[i][k] * L[j][k] ) / L[j][j]     for i > j, k = 0 .. j-1

Where:

- `A` is the symmetric positive-definite input matrix (`n x n`)
- `L` is the lower-triangular output, satisfying `A = L * L^T`

Each column depends on the columns computed before it, so the diagonal square roots and the off-diagonal updates are evaluated in order.

---

## 2. Why This Kernel Matters

Cholesky decomposition is one of the most widely used operations in numerical computing:

- **Solving linear systems**, where a positive-definite system is factored once and then solved cheaply for many right-hand sides
- **Least-squares and regression**, where the normal equations produce a positive-definite matrix to factor
- **State estimation and control**, where Kalman filters and related methods factor covariance matrices
- **Statistics and simulation**, where sampling correlated variables uses the Cholesky factor of a covariance matrix

Because it mixes reductions, a square root, and a division with strong dependencies between steps, it is a good measure of how well an architecture pipelines a structured, dependency-carrying computation.

---

## 3. Why EFF Hardware Performs Well

The Electron E1x runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

This fits fixed-point Cholesky decomposition well:

- The multiply-accumulate reductions that build each entry run **in parallel** across the Fabric
- Within a column, the off-diagonal entries are independent and are processed **several at a time**
- Partial sums and the current column stay **resident on the Fabric** instead of spilling to memory
- Fixed-point arithmetic keeps the datapath **simple**, so more of the Fabric is spent on useful work

The result is high arithmetic throughput at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## 4. Configurable Parameters

These definitions in `main.c` and `cholesky.h` control the benchmark. Change them to resize the problem or re-run it.

| Definition         | Default | Effect                                                                                                                                                                                                                                   |
| ------------------ | ------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `NUM_ITERATIONS`   | `1`     | How many times the kernel runs. Increase it to average out measurement noise when benchmarking.                                                                                                                                          |
| `INPUT_SIZE`       | `32`    | The matrix dimension `n` (the input is `n x n`). This sets the problem size: a larger value means a larger matrix and more work per run. The input is generated in `main.c` and made symmetric positive-definite before the kernel runs. |
| `FIXED_POINT_BITS` | `10`    | The number of fractional bits in the fixed-point format, defined in `cholesky.h`. It sets the trade-off between numeric range and precision. Changing it changes the scaling used throughout the kernel.                                 |
