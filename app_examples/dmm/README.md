# Dense Matrix Multiply (DMM)

This example computes a **dense matrix-matrix product** (`Z = A * B`) on the Electron E1 general-purpose processor. It is a compact, compute-dense linear-algebra kernel used to show how regular arithmetic maps onto the Fabric architecture.

---

## 1. Overview

### What is a Dense Matrix Multiply?
A dense matrix multiply takes an `n x m` matrix `A` and an `m x o` matrix `B` and produces an `n x o` matrix `Z`. Every element of the output is the dot product of one row of `A` with one column of `B`. "Dense" means every element is stored and processed (there is no sparsity to skip over).

### Mathematical Definition

    Z[i][j] = Σ A[i][k] * B[k][j]     for k = 0 .. m-1

Where:
- `A` is the left input matrix (`n` rows, `m` columns, row-major)
- `B` is the right input matrix (`m` rows, `o` columns, row-major)
- `Z` is the output matrix (`n` rows, `o` columns)

Each output element `Z[i][j]` costs `m` multiply-accumulate (MAC) operations, for `n * m * o` MACs total. For a square problem of size `n` this is O(n^3) work, which makes matrix multiply one of the most arithmetic-heavy kernels in common use.

---

## 2. Why This Kernel Matters

Matrix multiply is one of the most common operations in real systems:

- **Neural network inference**, where fully connected, convolutional, and attention layers are expressed as matrix-matrix products
- **Signal and image processing**, where transforms and filter banks apply a matrix to blocks of data
- **Scientific computing and solvers**, where linear systems and factorizations are built on repeated matrix products
- **Graphics and robotics**, where chained coordinate and pose transforms compose as matrix products

Because it is arithmetic-bound and highly regular, matrix multiply is a good measure of how well an architecture sustains dense computation while keeping data movement low.

---

## 3. Why EFF Hardware Performs Well

The Electron E1 runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

This fits dense matrix multiply well:

- Many multiply-accumulate operations run **in parallel** across the Fabric, so a whole tile of output elements advances at once
- The partial-sum accumulators for a tile stay **resident on the Fabric** instead of spilling to memory between steps
- Each loaded value from `A` and `B` is **reused across several accumulators** in the tile, so computation overlaps with loads and memory traffic stays low
- The regular, predictable access pattern maps cleanly onto a static dataflow graph, so very little energy goes to control overhead

The result is high arithmetic throughput at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## 4. Configurable Parameters

These definitions in `main.c` (and the shared matrix header it includes) control the benchmark. Change them to resize the problem or re-run it.

| Definition | Default | Effect |
|---|---|---|
| `NUM_ITERATIONS` | `1` | How many times the kernel runs. Increase it to average out measurement noise when benchmarking. |
| `MAT_REF_SIZE` | `32` | The matrix dimension (used for all of `n`, `m`, and `o`), defined in the shared `mat.h`. This sets the problem size: a larger value means larger matrices and more MACs per run, growing as O(n^3). The input matrices `mat_a` and `mat_b` are supplied from `mat.h`. |
