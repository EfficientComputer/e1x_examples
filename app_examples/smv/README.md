# Sparse Matrix-Vector Multiply (SMV)

This example computes a **sparse matrix-vector product** (`z = A * x`) on the Electron E1 general-purpose processor. The matrix is stored in compressed sparse row (CSR) format, so only the nonzero entries are processed. It is a good example of how irregular, indirect memory access maps onto the Fabric architecture.

---

## 1. Overview

### What is a Sparse Matrix-Vector Multiply?
A sparse matrix-vector multiply takes an `n x m` matrix `A` in which most entries are zero and multiplies it by a dense `m`-element vector `x` to produce an `n`-element vector `z`. Each output element is the dot product of one row of `A` with `x`, but only the nonzero entries of that row contribute, so the zeros are never stored or multiplied.

The matrix here is stored in CSR (Compressed Sparse Row) format. CSR keeps three arrays: a data array of the nonzero values, an indices array giving the column of each nonzero value, and a row-pointer array whose entry `i` marks where row `i` begins in the other two arrays. Row `i` therefore runs from `indptr[i]` up to `indptr[i+1]`.

### Mathematical Definition

    z[i] = Σ A[i][k] * x[k]     over the nonzero columns k in row i

Where:
- `A` is the input matrix (`n` rows, stored in CSR as data, indices, and row pointers)
- `x` is the dense input vector (`m` elements)
- `z` is the output vector (`n` elements)

Only the nonzero entries of each row take part, so the total work is proportional to the number of nonzeros, not to `n * m`.

---

## 2. Why This Kernel Matters

Sparse matrix-vector multiply is one of the most common operations wherever data is large but mostly empty:

- **Graph analytics**, where an adjacency matrix is sparse and repeated multiplies implement traversals and ranking algorithms
- **Scientific computing and solvers**, where iterative methods repeatedly apply a sparse system matrix to a vector
- **Recommendation and search**, where user-item and feature matrices are highly sparse
- **Machine learning**, where pruned or sparse layers skip the zero weights to save work

Because it is bound by irregular memory access rather than raw arithmetic, it is a good measure of how well an architecture handles indirect, data-dependent access patterns.

---

## 3. Why EFF Hardware Performs Well

The Electron E1 runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

Sparse matrix-vector multiply is dominated by irregular, indirect memory access and data-dependent control flow, which is exactly where a dataflow mapping helps:

- The index chasing that reads column positions **overlaps** with the multiply-accumulate arithmetic, hiding the latency of indirect loads
- Multiple output rows are computed **in parallel** across the Fabric, since each row is independent
- Row accumulators stay **resident on the Fabric** instead of spilling to memory
- Because only nonzero entries flow through the graph, **no work is spent** on the zeros

The result is high throughput on irregular data at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## 4. Configurable Parameters

These definitions in `main.c` (and the shared matrix header it includes) control the benchmark. Change them to resize the problem or re-run it.

| Definition | Default | Effect |
|---|---|---|
| `NUM_ITERATIONS` | `1` | How many times the kernel runs. Increase it to average out measurement noise when benchmarking. |
| `MAT_REF_SIZE` | `32` | The matrix dimension (number of rows and columns), defined in the shared `mat.h`. This sets the problem size. The sparse matrix `mat_a_sparse_*` and the dense `vector` are supplied from `mat.h`. |
