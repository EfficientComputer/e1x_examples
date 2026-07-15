# Recursive Factorial (Factorial)

This example computes a **factorial** using recursion on the Electron E1 general-purpose processor. It multiplies down through the integers from `n` to 1, showing how a recursive kernel maps onto the Fabric architecture.

---

## 1. Overview

### What is a factorial?
The factorial of a non-negative integer `n`, written `n!`, is the product of all the positive integers from 1 up to `n`. For example, `10!` is `10 * 9 * 8 * ... * 1`, which equals 3628800. Factorials grow very quickly and appear throughout counting problems, probability, and series expansions.

This implementation is recursive: it multiplies the running result by `n`, then calls itself on `n - 1`, stopping when `n` reaches 0. The kernel therefore exercises recursion and a chain of dependent multiplications.

### Mathematical Definition

    n! = n * (n-1) * (n-2) * ... * 1
    0! = 1

Expressed as the recurrence the kernel follows:

    factorial(n) = n * factorial(n - 1)
    factorial(0) = 1

Where each recursive call reduces `n` by one and multiplies its value into the running product until the base case is reached.

---

## 2. Why This Kernel Matters

Factorials and the recursion pattern behind them show up across many settings:

- **Combinatorics**, where counts of permutations and combinations are built from factorials
- **Probability and statistics**, where distributions use factorial terms
- **Series and numerical methods**, where factorial denominators appear in expansions
- **Teaching and benchmarking**, where recursion is a canonical example of a self-referential computation

As a benchmark it is a compact way to exercise recursion and a chain of dependent integer multiplications.

---

## 3. Why EFF Hardware Performs Well

The Electron E1 runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

The factorial is a sequence of dependent multiplications driven by recursion, and that pattern maps well onto the Fabric:

- The multiply-and-recurse steps are laid out as a **pipeline**, so the running product flows through without per-step instruction handling
- The running product stays **resident** near the compute elements instead of moving through memory
- The recursion follows a simple, regular structure that maps onto a compact dataflow graph, so very little energy goes to control overhead
- Advancing the counter and performing the multiplication **overlap**, keeping the compute elements busy

The result is efficient evaluation of a recursive integer computation at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## 4. Configurable Parameters

These definitions in `main.c` control the benchmark. Change them to re-run it.

| Definition | Default | Effect |
|---|---|---|
| `NUM_ITERATIONS` | `1` | How many times the kernel runs. Increase it to average out measurement noise when benchmarking. |
| `EXPECTED_OUTPUT` | `3628800` | The expected result the output is compared against for the pass/fail print. This is the value of `10!`. If you change the input, this must be updated to match the new correct result. |

The factorial argument is a hardcoded literal in `main.c` (the variable `n` is set to `10`). You can edit this literal to compute a different factorial, but if you do, update `EXPECTED_OUTPUT` to match the new correct result. Note that factorials grow quickly, so large arguments will overflow the 32-bit integer result.
