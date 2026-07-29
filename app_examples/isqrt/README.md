# Integer Square Root (isqrt)

This example computes the **integer square root** of a number on the Electron E1x general-purpose processor. It uses a digit-by-digit binary method that relies only on shifts, additions, and comparisons, showing how a compact integer kernel maps onto the Fabric architecture.

---

## Overview

### What is an integer square root?

The integer square root of a non-negative integer `x` is the largest integer `r` such that `r * r` is less than or equal to `x`. For example, the integer square root of 20000 is 141, because 141 squared is 19881 (at or below 20000) while 142 squared is 20164 (above it). The result is the exact square root rounded down to a whole number.

This implementation avoids both division and floating point. It resolves the answer two bits at a time, from the most significant bit downward, using only shifts, subtractions, additions, and comparisons.

### Mathematical Definition

The method builds the result one bit at a time. It first finds the highest power of four that is at or below `x`, then works downward. At each step it tests whether the current candidate bit can be included:

    t = x - r - q

Where:

- `q` is the current power-of-four place value, halved by two bit positions each step
- `r` is the partial result accumulated so far
- if `t` is at or above zero, the bit is accepted: `x` is reduced to `t` and `q` is added into `r`
- `r` is shifted right by one bit each step so it stays aligned with the shrinking place value

When `q` reaches its final step, `r` holds the exact integer square root.

---

## Why This Kernel Matters

Integer square roots are needed wherever a magnitude must be computed without floating-point hardware:

- **Fixed-point signal processing**, where amplitudes and magnitudes are computed on integer data
- **Distance and geometry**, where Euclidean lengths are derived from sums of squares
- **Embedded and sensor firmware**, where small controllers lack a floating-point unit
- **Graphics and physics**, where normalization and radius checks use square roots
- **Statistics on integer data**, where standard deviations and norms are computed

It is a good example of a small, branch-driven integer kernel built entirely from shifts, adds, and comparisons.

---

## Why Efficient Hardware Performs Well

The Electron E1x runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

The digit-by-digit method is a short sequence of dependent integer steps, and that pattern maps well onto the Fabric:

- Each bit-resolving step is laid out as a **pipeline** stage, so values flow through without per-step instruction handling
- The partial result and remaining value stay **resident** near the compute elements instead of moving through memory
- The subtract-compare-and-select logic that decides each bit flows directly between compute elements as a small dataflow graph
- Using only shifts, adds, and comparisons keeps the operation simple, so very little energy goes to control overhead

The result is an exact square root with no division or floating point at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## Configurable Parameters

This example has no problem-size macros. The kernel computes the integer square root of whatever value it is given.

| Definition       | Default | Effect                                                                                                                 |
| ---------------- | ------- | ---------------------------------------------------------------------------------------------------------------------- |
| `NUM_ITERATIONS` | `1`     | This is how many times the first test case runs (in `main.c`). Increase it to average out measurement noise when benchmarking. |

The inputs are hardcoded literals in `main.c`: the kernel is called on `16`, `1764`, and `20000`, and each result is checked against its expected value (`4`, `42`, and `141`). You can edit these literals to test different values, but if you do, update the expected results in the same checks to match the new correct answers.
