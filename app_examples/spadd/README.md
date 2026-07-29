# Sparse Matrix Addition (SpAdd)

This example adds two **sparse matrices** (`Z = A + B`) on the Electron E1x general-purpose processor, producing a dense output matrix. Both inputs are stored in compressed sparse row (CSR) format, so only the nonzero entries are processed. It shows how irregular, data-dependent access maps onto the Fabric architecture.

---

## Overview

### What is a Sparse Matrix Addition?

A sparse matrix addition takes two `n x n` matrices `A` and `B`, both mostly zero, and adds them element by element to produce a matrix `Z`. Where a position is nonzero in only one input, that value carries through; where it is nonzero in both, the two values are summed. The result is written out in dense form.

Both inputs are stored in CSR (Compressed Sparse Row) format. CSR keeps three arrays: a data array of the nonzero values, a coordinate array giving the column of each nonzero value, and a row-pointer array whose entry `i` marks where row `i` begins in the other two arrays. Adding the two matrices amounts to merging the nonzero entries of each row of `A` with the matching row of `B`.

### Mathematical Definition

    Z[i][j] = A[i][j] + B[i][j]     for every position (i, j)

Where:

- `A` and `B` are the input matrices (both `n x n`, stored in CSR as data, coordinates, and row pointers)
- `Z` is the dense output matrix (`n x n`)

Only positions that are nonzero in `A` or `B` contribute, so the work is proportional to the total number of nonzeros in the two inputs rather than to `n * n`.

---

## Why This Kernel Matters

Adding sparse matrices is a building block wherever large, mostly-empty operators are combined:

- **Graph analytics**, where combining sparse adjacency or weight matrices merges or updates graphs
- **Scientific computing and solvers**, where sparse system matrices are assembled and updated term by term
- **Finite-element and mesh methods**, where local contributions are accumulated into a global sparse operator
- **Machine learning**, where sparse gradient or weight updates are applied to sparse parameters

Because it is bound by irregular access and the data-dependent merging of two sparse streams, it is a good measure of how well an architecture handles indirect access and branchy control flow.

---

## Why Efficient Hardware Performs Well

The Electron E1x runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

This kernel is dominated by irregular, indirect memory access and data-dependent control flow, which is exactly where a dataflow mapping helps by overlapping index chasing with arithmetic:

- The merge that lines up nonzero column positions **overlaps** with the addition, hiding the latency of indirect loads
- Each output row is independent, so multiple rows are processed **in parallel** across the Fabric
- Because the row merge does not depend on the order in which stores land, the memory-ordering constraints can be **relaxed**, exposing more parallelism
- Because only nonzero entries flow through the graph, **no work is spent** on positions that are zero in both inputs

The result is high throughput on irregular data at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## Configurable Parameters

These definitions in `main.c` (and the shared matrix header it includes) control the benchmark. Change them to resize the problem or re-run it.

| Definition       | Default | Effect                                                                                                                                                                                                                |
| ---------------- | ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `NUM_ITERATIONS` | `1`     | This is how many times the kernel runs. Increase it to average out measurement noise when benchmarking.                                                                                                                       |
| `MAT_REF_SIZE`   | `32`    | This is the matrix dimension (number of rows and columns). This sets the problem size and the size of the dense output `mat_z`. The sparse inputs `mat_a_sparse_*` and `mat_b_sparse_*` are supplied from the shared `mat.h`. |
