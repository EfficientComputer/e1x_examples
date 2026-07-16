# Sparse Vector Operations (SVec)

This example performs **operations on sparse vectors** on the Electron E1x general-purpose processor. Each vector is stored as coordinate-value pairs, listing only its nonzero entries. The core operation is a stream-join multiply of two sparse vectors, which is also the basis of a sparse dot product. It shows how irregular, data-dependent access maps onto the Fabric architecture.

---

## 1. Overview

### What is a Sparse Vector Operation?

A sparse vector stores only its nonzero entries, each as a coordinate (its position) and a value. Two sparse vectors are combined by walking both coordinate lists together and acting only where their coordinates match. The stream-join multiply produces a new sparse vector whose entries are the products at the shared coordinates. Summing those products gives the sparse dot product of the two vectors.

The example supports two ways of holding the vectors. A fixed representation uses statically sized arrays with hand-picked contents so the result can be checked against a known answer. A dynamic representation allocates the vectors and fills them with random nonzeros. The default build uses the fixed representation.

### Mathematical Definition

For the stream-join multiply of sparse vectors `a` and `b`, the result `c` has an entry at every coordinate `k` where both inputs are nonzero:

    c[k] = a[k] * b[k]     for every coordinate k present in both a and b

The sparse dot product then reduces those products to a single scalar:

    dot(a, b) = Σ a[k] * b[k]     over coordinates k present in both a and b

Only shared coordinates contribute, so the work depends on the overlap of the two vectors rather than on their full length.

---

## 2. Why This Kernel Matters

Sparse vector operations are a core building block wherever data is high-dimensional but mostly empty:

- **Search and information retrieval**, where sparse term vectors are scored against one another
- **Recommendation systems**, where sparse user and item profiles are compared by dot product
- **Machine learning**, where sparse feature vectors and sparse activations are combined
- **Graph analytics**, where sparse frontier and neighbor vectors are intersected during traversal

Because it is bound by irregular access and the data-dependent matching of two sparse streams, it is a good measure of how well an architecture handles indirect access and branchy control flow.

---

## 3. Why EFF Hardware Performs Well

The Electron E1x runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

This kernel is dominated by irregular, indirect memory access and data-dependent control flow, which is exactly where a dataflow mapping helps by overlapping index chasing with arithmetic:

- The stream-join that advances the two coordinate lists **overlaps** with the multiply, hiding the latency of indirect loads
- The running product and reduction stay **resident on the Fabric** instead of spilling to memory
- The comparisons that decide which pointer advances become part of the dataflow graph, so **little energy** goes to control overhead
- Because only the matching coordinates flow through the graph, **no work is spent** on positions that are zero in either vector

The result is high throughput on irregular data at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## 4. Configurable Parameters

These definitions in `svec.h` and `main.c` control the benchmark. Change them to resize the problem or switch the vector representation.

| Definition       | Default   | Effect                                                                                                                                                                                                                                         |
| ---------------- | --------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `SVEC_FIXED_VEC` | defined   | Uses the fixed, hand-picked input vectors so the result can be checked against a known answer.                                                                                                                                                 |
| `SVEC_DYNAMIC`   | undefined | When defined, allocates the vectors dynamically instead of using the static arrays. Leave undefined to use the static representation.                                                                                                          |
| `SVEC_MAX_NNZS`  | `100`     | The maximum number of nonzero entries a vector can hold. Sets the size of the static coordinate and value arrays.                                                                                                                              |
| `MAX_VEC_NNZS`   | `100`     | The nonzero cap used when filling a random vector. Set equal to `SVEC_MAX_NNZS`.                                                                                                                                                               |
| `MAX_COORD`      | `100`     | The largest coordinate a random vector can use, that is, the effective length of the vector space.                                                                                                                                             |
| `NNZ_LIKELIHOOD` | `8`       | For random vectors, the probability that any given coordinate is nonzero is 1 / `NNZ_LIKELIHOOD`. A larger value makes the vectors sparser.                                                                                                    |
| `V1LEN`          | `7`       | Number of nonzeros in the first fixed input vector.                                                                                                                                                                                            |
| `V2LEN`          | `11`      | Number of nonzeros in the second fixed input vector.                                                                                                                                                                                           |
| `RESLEN`         | `2`       | Number of nonzeros in the expected result. With `SVEC_FIXED_VEC`, the output is checked against `res_index` and `res_val`. If you change the fixed inputs, these expected values and `RESLEN` must be updated to match the new correct result. |
