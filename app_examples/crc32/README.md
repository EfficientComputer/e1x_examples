# CRC-32 Checksum (CRC-32)

This example computes a **CRC-32 checksum** over a block of data on the Electron E1x general-purpose processor. It processes the message one byte at a time to produce a 32-bit checksum, showing how a bitwise, sequential kernel maps onto the Fabric architecture.

---

## Overview

### What is CRC-32?

A cyclic redundancy check (CRC) is a checksum used to detect accidental changes to data. CRC-32 produces a 32-bit value by treating the message as a large binary number and computing its remainder against a fixed generator polynomial using carry-less (XOR-based) arithmetic. The same input always yields the same checksum, so a receiver can recompute the CRC and compare it to detect corruption.

This implementation uses the standard reflected CRC-32 with polynomial `0xEDB88320`. It starts from an all-ones value, folds in each message byte, and inverts the result at the end.

### Mathematical Definition

For each byte of the message, the running CRC value is updated bit by bit. Eight times per byte:

    mask = -(crc & 1)
    crc  = (crc >> 1) ^ (0xEDB88320 & mask)

Where:

- `crc` starts at `0xFFFFFFFF` and each incoming byte is XORed into its low bits before the eight shift steps
- `0xEDB88320` is the reflected form of the standard CRC-32 generator polynomial
- `mask` is all ones when the low bit is set and all zeros otherwise, so the polynomial is applied only when needed
- the final checksum is the bitwise complement (`~crc`) of the running value after the last byte

---

## Why This Kernel Matters

CRC-32 is one of the most widely deployed integrity checks in computing:

- **Network protocols**, where frames and packets carry a CRC so corruption in transit is caught
- **Storage and file formats**, where archives and disk sectors are checked for silent errors
- **Firmware and boot images**, where an image is validated before it is trusted
- **Embedded communication buses**, where short messages are guarded against noise
- **Data pipelines**, where blocks are verified as they move between stages

It is a good example of a bitwise, control-driven kernel: the work is shifts, masks, and XORs rather than multiply-heavy arithmetic.

---

## Why Efficient Hardware Performs Well

The Electron E1x runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

CRC-32 is a tight loop of bit operations over a stream of bytes, and that pattern maps well onto the Fabric:

- The eight shift-and-XOR steps per byte are laid out as a **pipeline**, so bits flow through without per-step instruction handling
- The running CRC value stays **resident** near the compute elements instead of spilling to memory
- Loading the next byte **overlaps** with folding the current byte, hiding memory latency behind useful work
- The mask-and-XOR logic is a small, regular dataflow graph, so very little energy goes to control overhead

The result is efficient integrity checking at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## Configurable Parameters

These definitions in `main.c` control the benchmark. Change them to resize the problem or re-run it.

| Definition       | Default | Effect                                                                                                                                                                                                          |
| ---------------- | ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `NUM_ITERATIONS` | `1`     | This is how many times the kernel runs. Increase it to average out measurement noise when benchmarking.                                                                                                                 |
| `LENGTH`         | `512`   | This is the size in bytes of the message buffer. This sets the problem size: a larger value means more bytes to fold in per run. The last byte is kept as a null terminator, so the checksum covers `LENGTH - 1` bytes. |

The message is filled with pseudo-random data generated in `main.c` from a fixed seed (`123456789`), so the input is the same on every run. Changing `LENGTH` or the seed changes the input data and therefore the result.
