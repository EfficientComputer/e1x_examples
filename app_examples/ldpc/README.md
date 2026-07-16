# Low-Density Parity-Check Code (LDPC)

This example encodes and decodes a Low-Density Parity-Check error-correcting code on the Electron E1x general-purpose processor. It encodes a message, injects bit errors, and iteratively decodes it back to the original codeword, a workload that is central to modern communications and storage.

---

## 1. Overview

### What is an LDPC Code?

An LDPC code is a linear error-correcting code defined by a sparse parity-check matrix. Encoding multiplies the message bits by a generator matrix (modulo 2) to add redundant parity bits. Decoding uses iterative message passing between two sets of nodes (variable nodes, which represent codeword bits, and check nodes, which represent parity constraints) to recover the original data even when some bits are corrupted. This example uses the min-sum approximation of belief propagation, operating on log-likelihood ratios (LLRs) that express how confident the decoder is about each bit.

### Mathematical Definition

Encoding computes the codeword `c` from the message `m` using the generator matrix `G`, where all arithmetic is modulo 2 (that is, addition is exclusive-or):

    c[j] = XOR over i of ( G[i][j] * m[i] )

Decoding is iterative. On each iteration the decoder refines the LLR for every bit by exchanging messages along the sparse parity-check matrix `H`, using the min-sum rule at each check node:

    check-to-variable message ~= ( product of signs ) * ( minimum of magnitudes )

It stops early once the current bit decisions satisfy every parity check (the syndrome is all zero) or once the maximum iteration count is reached.

---

## 2. Why This Kernel Matters

LDPC codes deliver near-optimal error correction with practical decoding cost, so they are built into many communications and storage standards:

- 5G New Radio, where LDPC protects the data channel
- Wi-Fi (802.11n/ac/ax), the family this example's code is drawn from
- Satellite and deep-space links, where signals are weak and errors are common
- High-speed wired links such as 10GBASE-T Ethernet
- Solid-state drives and other storage, where LDPC extends reliability and endurance

Because decoding combines sparse, irregular memory access with iterative message passing, it is a strong measure of how well an architecture handles structured but non-uniform workloads.

---

## 3. Why EFF Hardware Performs Well

The Electron E1x runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

This fits LDPC decoding well:

- The parity-check and generator matrices are stored in compressed sparse row (CSR) form, so only the non-zero connections are traversed and no work is spent on empty entries
- The many independent check nodes are processed in parallel across the Fabric, since each one operates on its own small set of connected bits
- The forward and backward min-sum scans within a check node are decomposed into sign and magnitude parts, which keeps the recurrence tight and streaming
- Message and LLR values pass directly between compute elements from one processing step to the next instead of spilling to memory
- Early termination stops the iteration as soon as a valid codeword is found, avoiding wasted work

The result is high error-correction throughput at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## 4. Configurable Parameters

These definitions and variables control the benchmark. The problem-size constants are `#define`s in `main.c`; the decoder controls are local variables in the `test()` function. Change them to use a different code, adjust decoding effort, or re-run the benchmark.

| Definition                    | Default | Effect                                                                                                                                                        |
| ----------------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `NUM_ITERATIONS`              | `1`     | How many times the decode kernel runs. Increase it to average out measurement noise when benchmarking.                                                        |
| `MESSAGE_SIZE`                | `336`   | The number of message (information) bits. Set by the 802.11 rate-1/2 code used here.                                                                          |
| `CODEWORD_SIZE`               | `672`   | The number of codeword bits (message plus parity). Set by the same rate-1/2 code.                                                                             |
| `PARITY_NNZ`                  | `2184`  | The number of non-zero entries in the parity-check matrix H, in CSR form.                                                                                     |
| `GENERATOR_NNZ`               | `8736`  | The number of non-zero entries in the generator matrix G, in CSR form.                                                                                        |
| `MAX_CHECK_NODE_NNZ`          | `8`     | The maximum number of bits connected to any single check node, a property of this code's structure.                                                           |
| `decoder_max_iterations`      | `16`    | The maximum number of decoding iterations before the decoder gives up. Higher values allow correcting more difficult error patterns at the cost of more work. |
| `decoder_early_term_possible` | `true`  | Whether the decoder stops as soon as a valid codeword is found. Disabling it forces every iteration to run.                                                   |

The parity-check and generator matrices (`H_80211`, `G_80211`, and their CSR views) are fixed data for the 802.11 rate-1/2 code and are generated to match these size constants. The test message and the injected bit errors are set in `main.c`. If you change the code or its dimensions, the matrix data and the size constants must be updated together to match the new correct result.
