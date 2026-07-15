# Dense Matrix-Vector Multiply (DMV)

This example computes a **dense matrix-vector product** (`z = A · b`) on the Electron E1 general-purpose processor. It is a compact, compute-dense linear-algebra kernel used to show how regular arithmetic maps onto the Fabric architecture.

---

## 1. Overview

### What is a Dense Matrix-Vector Multiply?
A dense matrix-vector multiply takes an `n × m` matrix `A` and an `m`-element vector `b` and produces an `n`-element vector `z`. Every element of the output is the dot product of one row of `A` with the vector `b`. "Dense" means every element is stored and processed (there is no sparsity to skip over).

### Mathematical Definition

    z[i] = Σ A[i][j] * b[j]     for j = 0 .. m-1

Where:
- `A` is the input matrix (`n` rows, `m` columns, row-major)
- `b` is the input vector (`m` elements)
- `z` is the output vector (`n` elements)

Each output element `z[i]` costs `m` multiply-accumulate (MAC) operations, for `n × m` MACs total.

---

## 2. Why This Kernel Matters

Matrix-vector multiply is one of the most common operations in real systems:

- **Neural network inference**, where fully connected and attention-projection layers are matrix-vector (or matrix-matrix) products
- **Signal and control systems**, where state updates and filters apply a matrix to a measurement vector
- **Scientific computing and solvers**, where iterative methods repeatedly apply a matrix to a vector
- **Graphics and robotics**, where coordinate and pose transforms are small matrix-vector products

Because it is memory-bound at large sizes and arithmetic-bound at small sizes, it is a good measure of how well an architecture balances computation with data movement.

---

## 3. Why EFF Hardware Performs Well

The Electron E1 runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

This fits dense matrix-vector multiply well:

- Many multiply-accumulate operations run **in parallel** across the Fabric
- Row accumulators stay **resident on the Fabric** instead of spilling to memory
- Computation **overlaps with memory loads**, hiding load latency
- The regular, predictable access pattern maps cleanly onto a static dataflow graph, so very little energy goes to control overhead

The result is high arithmetic throughput at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## 4. Configurable Parameters

These definitions in `main.c` (and the shared matrix header it includes) control the benchmark. Change them to resize the problem or re-run it.

| Definition | Default | Effect |
|---|---|---|
| `NUM_ITERATIONS` | `1` | How many times the kernel runs. Increase it to average out measurement noise when benchmarking. |
| `MAT_REF_SIZE` | `32` | The matrix dimension (both `n` and `m`), defined in the shared `mat.h`. This sets the problem size: a larger value means a larger matrix and vector and more MACs per run. The input matrix `mat_a` and `vector` are supplied from `mat.h`. |
