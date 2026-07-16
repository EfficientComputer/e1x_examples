# Newton-Raphson Square Root (sqrt_newton)

This example computes **square roots** using Newton-Raphson iteration on the Electron E1x general-purpose processor. It refines a floating-point estimate through a fixed number of iterations, showing how an iterative numerical kernel maps onto the Fabric architecture.

---

## 1. Overview

### What is Newton-Raphson square root?

Newton-Raphson is a general method for finding where a function equals zero by repeatedly improving a guess. To compute the square root of `a`, the method solves `f(x) = x * x - a = 0`. Starting from an initial estimate, each iteration moves the estimate closer to the true root, and the number of correct digits roughly doubles per step, so a handful of iterations is enough for good accuracy.

This implementation runs a fixed number of iterations on floating-point values and returns the refined estimate. It handles the special cases of zero and negative input directly.

### Mathematical Definition

Each iteration updates the current estimate `x` using the Newton-Raphson step for `f(x) = x * x - a`:

    x_next = x - (x * x - a) / (2 * x)

Where:

- `a` is the value whose square root is being computed
- `x` is the current estimate, refined toward the square root of `a`
- the denominator `2 * x` is the derivative of `x * x - a`
- the loop runs a fixed number of times, after which `x` is returned as the result

---

## 2. Why This Kernel Matters

Square roots computed by iteration appear across numerical and signal-processing code:

- **Vector normalization**, where lengths are needed to scale directions in graphics and physics
- **Signal processing**, where magnitudes are taken from sums of squares
- **Statistics**, where standard deviations and norms are computed
- **Control and estimation**, where gains and covariances involve square roots
- **Scientific computing**, where roots appear inside larger iterative solvers

It is a good example of an iterative floating-point kernel where the same small update is applied repeatedly, and accuracy is traded against the number of iterations.

---

## 3. Why EFF Hardware Performs Well

The Electron E1x runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

Newton-Raphson is a short chain of arithmetic repeated a fixed number of times, and that pattern maps well onto the Fabric:

- The multiply, subtract, and divide of each iteration are laid out as a **pipeline**, so values flow through without per-step instruction handling
- The running estimate stays **resident** near the compute elements instead of spilling to memory
- Because the iteration count is fixed, the whole update chain maps onto a **static dataflow graph**, so very little energy goes to control overhead
- Successive square-root computations can be **overlapped**, keeping the compute elements busy

The result is fast convergence to an accurate estimate at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## 4. Configurable Parameters

These definitions in `sqrt_newton.c` (which holds both the kernel and the test driver) control the benchmark. Change them to resize the problem or re-run it.

| Definition | Default | Effect                                                                                                                                                   |
| ---------- | ------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `SQRT_QTY` | `10`    | How many values have their square root computed. The driver computes the square root of each integer from 0 up to `SQRT_QTY - 1` and checks each result. |

The number of refinement iterations is a hardcoded literal in the `sqrt_newton` function (the loop runs `5` times). Increasing it improves accuracy at the cost of more work; decreasing it does the reverse. The expected results are stored in the `actual_sqrt` array in `sqrt_newton.c`, compared with a tolerance in `float_equality`. If you change `SQRT_QTY` or the set of inputs, the `actual_sqrt` array must be updated to match the new correct values.
