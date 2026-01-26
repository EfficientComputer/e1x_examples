# Optimized LDPC Encoder/Decoder

This example provides an **optimized LDPC (Low-Density Parity-Check) encoder and decoder** tailored for Efficient Computer (EFF) hardware.  
It is designed to showcase **real-world error correction performance** and highlight the advantages of EFF's architecture for communications and storage workloads.

---

## 1. Overview

### What is LDPC?
**LDPC (Low-Density Parity-Check)** codes are a class of linear error-correcting codes widely used in modern communications systems. They provide near-Shannon-limit performance with practical decoding complexity.

Common use cases include:
- 5G NR (New Radio) wireless communications
- WiFi (802.11n/ac/ax)
- DVB-S2 satellite broadcasting
- 10GBASE-T Ethernet
- Solid-state drive (SSD) error correction
- Deep-space communications

---

### How LDPC Works

**Encoding**:
- Multiply message bits by a sparse generator matrix (mod 2)
- Produces codeword with redundant parity bits

**Decoding** (Belief Propagation / Min-Sum):
- Initialize log-likelihood ratios (LLRs) from received signal
- Iteratively pass messages between variable nodes and check nodes
- Converge to error-corrected codeword

The decoder uses the **min-sum approximation** for efficient hardware implementation.

---

## 2. Why This Kernel Matters

LDPC decoding is a **critical and challenging benchmark** because it:

- Involves **sparse matrix operations** with irregular memory access
- Requires **iterative message passing** with loop-carried dependencies
- Has **variable convergence** based on channel conditions
- Exhibits **fine-grained parallelism** across nodes but with dependencies between iterations
- Is **latency-sensitive** in real-time communications

These characteristics make it an excellent test of **architectural efficiency** for communications workloads.

---

## 3. Why EFF Hardware Performs Well

This implementation is optimized to:

- Use **CSR (Compressed Sparse Row)** format for efficient sparse matrix traversal
- Decompose **min-sum computation** into parallel forward/backward scans
- Exploit EFF's ability to **parallelize across check nodes** while preserving iteration order
- Keep message arrays and accumulators in **registers** where possible
- Use **early termination** when valid codeword is detected
- Minimize carry dependencies with **explicit dependency decomposition**

Key optimization techniques:
- Separate forward (F) and backward (B) accumulators per check node enable parallel processing
- Sign and magnitude computed independently to reduce initiation interval
- Memory ordering pragmas allow aggressive instruction scheduling

---

## 4. Code Structure

### Core Files

- **`ldpc.c`**
  - LDPC encoder and decoder implementation
  - `encode(...)`: Sparse matrix multiply for encoding
  - `decode(...)`: Min-sum belief propagation decoder
  - `initialize_check_node_estimates(...)`: Initialize variable-to-check messages
  - `compute_f_and_b(...)`: Forward/backward min-sum scans (parallelized)
  - `update_c2v(...)`: Check-to-variable message updates
  - `update_message_nodes0/1(...)`: Variable node updates
  - `is_codeword(...)`: Syndrome check for early termination

- **`ldpc.h`**
  - Function declarations for encode/decode APIs

- **`ldpc_utils.h`**
  - Data structure definitions
  - `csr_matrix_t`: Compressed sparse row matrix format
  - `sparse_matrix_t`: Combined row/column views
  - `decode_tmp_allocations_t`: Working memory for decoder

- **`parity_80211.c`**
  - 802.11 WiFi rate-1/2 parity check matrix (336×672)
  - Pre-computed CSR format for efficient traversal

- **`generator_80211.c`**
  - 802.11 WiFi generator matrix
  - Used for encoding test messages

- **`main.c`**
  - Application entry point
  - Creates test message and encodes it
  - Introduces bit errors (erasures and flips)
  - Runs decoder and validates correction
  - Reports iteration count and PASS/FAIL status

- **`CMakeLists.txt`**
  - Build configuration
  - Defines the `ldpc` app and its build targets (sim / fabric / scalar)

---

## 5. Specification

- **Code**: (336, 672) rate-1/2 LDPC from 802.11 WiFi standard
- **Message Size**: 336 bits
- **Codeword Size**: 672 bits
- **Decoder Algorithm**: Min-sum approximation for belief propagation
- **Max Iterations**: 16 (configurable)
- **Early Termination**: Enabled (stops when valid codeword detected)
- **Input Format**: int32 log-likelihood ratios (LLR)
  - Negative LLR indicates bit is likely 1
  - Positive LLR indicates bit is likely 0

---

## 6. Build Instructions

### Prerequisites
- CMake
- EFF SDK environment / toolchain
- Math library (linked automatically)

### Build Steps (typical EFF SDK flow)
Refer to the main build instructions.

The build compiles:
- LDPC encoder and decoder
- 802.11 WiFi code matrices
- Error injection and validation in `main.c`

---

## 7. How to Run and Test

### Test Scenario
The test program:
1. Creates a test message with specific bit patterns
2. Encodes the message using the generator matrix
3. Converts to LLRs (±1000 representing high confidence)
4. Introduces errors:
   - 3 bit erasures (LLR set to 0)
   - 10+ bit flips (LLR sign inverted)
5. Runs the decoder
6. Validates that all errors are corrected

### Running
Run the built binary using your normal EFF workflow (sim or fabric). The app:
1. Encodes test message
2. Injects bit errors
3. Runs min-sum decoder (profiled region)
4. Compares decoded output to original codeword
5. Reports iteration count and PASS/FAIL status

### Correctness Checking
Correctness is validated by comparing each decoded bit:

```
_test_codeword_decision[i] == _test_codeword[i]
```

for all 672 bits in the codeword.

### Expected Output
```
LDPC error corrected in 2 iterations
[ldpc] PASS
```

---

## 8. Performance Benchmarking

This example is intended to highlight:
- Throughput for **iterative belief propagation**
- Efficiency of **sparse matrix operations**
- Performance with **irregular memory access patterns**
- EFF's advantages for **communications workloads**

Performance can be measured using:
- Cycle counters
- EFF simulator statistics
- Hardware profiling tools
- Iterations to convergence

---

## 9. Why This Example Is Useful

LDPC decoding is representative of:
- 5G NR physical layer processing
- WiFi/Bluetooth error correction
- Satellite communications
- Storage system reliability

Efficient execution of LDPC decoding directly impacts:
- Wireless throughput and latency
- Power consumption in mobile devices
- Error correction capability under poor channel conditions
- Storage system reliability and endurance

---

## 10. Summary

- Demonstrates an optimized **LDPC encoder and decoder**
- Implements **802.11 WiFi rate-1/2 code** (336, 672)
- Uses **min-sum approximation** for efficient decoding
- Features **early termination** for faster convergence
- Validates correctness with bit error injection
- Highlights **communications workload efficiency** on EFF hardware
- Relevant to 5G, WiFi, satellite, and storage applications
