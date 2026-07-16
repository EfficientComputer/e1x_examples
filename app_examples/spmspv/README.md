# Sparse Matrix-Sparse Vector Multiply (SpMSpV)

This example computes the product of a **sparse matrix and a sparse vector** (`z = A * x`) on the Electron E1x general-purpose processor, producing a dense output vector. The matrix is stored in compressed sparse row (CSR) format and the vector is stored as index-value pairs, so only nonzero entries take part. It shows how irregular access on both operands maps onto the Fabric architecture.

---

## 1. Overview

### What is a Sparse Matrix-Sparse Vector Multiply?

A sparse matrix-sparse vector multiply takes an `n x m` matrix `A` in which most entries are zero and multiplies it by a vector `x` that is also mostly zero, producing an `n`-element output vector `z`. Each output element is the dot product of one row of `A` with `x`. Because both the row and the vector are sparse, the dot product is computed as an intersection: only the positions where both operands have a nonzero value contribute.

The matrix is stored in CSR (Compressed Sparse Row) format. CSR keeps three arrays: a data array of the nonzero values, an indices array giving the column of each nonzero value, and a row-pointer array whose entry `i` marks where row `i` begins in the other two arrays. The sparse vector is stored as a list of its nonzero indices and the matching values.

### Mathematical Definition

    z[i] = Σ A[i][k] * x[k]     over columns k where both A[i][k] and x[k] are nonzero

Where:

- `A` is the input matrix (`n` rows, stored in CSR as data, indices, and row pointers)
- `x` is the sparse input vector (stored as nonzero indices and values)
- `z` is the dense output vector (`n` elements)

Each row's dot product walks the row and the vector together and accumulates only where their indices match, so the work depends on the overlap of the two sparse operands.

---

## 2. Why This Kernel Matters

Sparse matrix times sparse vector shows up wherever both the operator and the data are large but mostly empty:

- **Graph analytics**, where multiplying a sparse adjacency matrix by a sparse frontier vector advances a breadth-first traversal
- **Recommendation and search**, where a sparse feature vector is scored against a sparse item or user matrix
- **Scientific computing**, where sparse systems are applied to sparse states or residuals
- **Machine learning**, where sparse activations are combined with sparse weight matrices

Because it is bound by irregular access and the data-dependent matching of two sparse streams, it is a good measure of how well an architecture handles indirect access and branchy control flow.

---

## 3. Why EFF Hardware Performs Well

The Electron E1x runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

This kernel is dominated by irregular, indirect memory access and data-dependent control flow, which is exactly where a dataflow mapping helps by overlapping index chasing with arithmetic:

- The stream-join that matches nonzero column indices **overlaps** with the multiply-accumulate arithmetic, hiding the latency of indirect loads
- Multiple output rows are processed **in parallel** across the Fabric, since each row is an independent dot product
- Row accumulators stay **resident on the Fabric** instead of spilling to memory
- Because only the intersecting nonzeros flow through the graph, **no work is spent** on positions that are zero in either operand

The result is high throughput on doubly-sparse data at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## 4. Configurable Parameters

These definitions in `main.c` (and the shared matrix header it includes) control the benchmark. Change them to resize the problem or re-run it.

| Definition       | Default | Effect                                                                                                                                                                                                                  |
| ---------------- | ------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `NUM_ITERATIONS` | `1`     | How many times the kernel runs. Increase it to average out measurement noise when benchmarking.                                                                                                                         |
| `MAT_REF_SIZE`   | `32`    | The matrix dimension (number of rows and columns), defined in the shared `mat.h`. This sets the problem size. The sparse matrix `mat_a_sparse_*` and the sparse vector `vector_sparse_csc_*` are supplied from `mat.h`. |
