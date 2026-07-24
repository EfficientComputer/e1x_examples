# Biquad IIR Filter (Biquad)

This example runs a biquad (second-order) IIR filter over a block of audio samples on the Electron E1 general-purpose processor. It is a compact, feedback-based digital signal processing kernel used to show how a filter with loop-carried dependencies maps onto the Fabric architecture.

---

## 1. Overview

### What is a Biquad Filter?
A biquad is a second-order recursive (infinite impulse response, or IIR) digital filter. It computes each output sample from the two most recent input samples and the two most recent output samples, using five coefficients. Larger filters are built by chaining several biquad sections in series. In this example the filter is configured as a low-pass filter and runs in Q15 fixed-point arithmetic.

### Mathematical Definition

    y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2]

Where:
- `x[n]` is the current input sample, and `x[n-1]`, `x[n-2]` are the two previous inputs
- `y[n]` is the current output sample, and `y[n-1]`, `y[n-2]` are the two previous outputs
- `b0`, `b1`, `b2` are the feed-forward (zero) coefficients
- `a1`, `a2` are the feedback (pole) coefficients

The feedback terms create a loop-carried dependence: each output depends on outputs already produced, so samples cannot simply be computed independently in any order.

---

## 2. Why This Kernel Matters

Biquad filters are a staple of real-time signal processing and show up wherever a signal needs to be shaped:

- Audio processing, such as equalizers, tone controls, and crossover networks
- Sensor and instrumentation front ends, where noise is removed before further analysis
- Communications, where signals are band-limited or shaped before transmission
- Control systems, where measurements are smoothed or compensated
- Always-on and wearable devices, where continuous filtering must run at very low energy

Because they are efficient, numerically well behaved, and easy to cascade, biquads are the default way to implement many practical filters.

---

## 3. Why EFF Hardware Performs Well

The Electron E1 runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

This fits the biquad filter well:

- The five multiply-accumulate operations per sample map onto a fixed dataflow pattern that runs the same way for every sample
- The filter state (the previous inputs and outputs) stays resident on the Fabric and passes directly from one sample to the next instead of spilling to memory
- The loop-carried feedback path is expressed explicitly in the dataflow graph, so the compiler can keep the recurrence tight and streaming
- Input samples flow in and output samples flow out continuously, overlapping computation with data movement
- The regular, predictable structure means very little energy goes to control overhead

The result is steady, low-latency filtering at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## 4. Configurable Parameters

These definitions and constants in `main.c` control the benchmark. Change them to resize the problem, retune the filter, or re-run it.

| Definition | Default | Effect |
|---|---|---|
| `NUM_ITERATIONS` | `1000` | How many times the kernel runs inside the AON4-marked measured region. The kernel is short, so it is looped to exceed the EVK sample period; the measured runtime and energy are divided by this count in post-processing (pass the same value to `process-energy.py --iters`). Increase it if the region does not run comfortably longer than the sample period. |
| `INPUT_SIZE` | `16384` | The number of samples in the input block. This sets the problem size. If you change it, the build must be updated so the sample input and expected output data are regenerated for the new size. |
| `zeros_coeffs` | `{0.00306..., 0.00612..., 0.00306...}` | The three feed-forward (zero) coefficients `b0`, `b1`, `b2`. They define the filter response, set here for a low-pass filter centered near 300 Hz with Q around 0.72. |
| `poles_coeffs` | `{-1.84008..., 0.85233...}` | The two feedback (pole) coefficients `a1`, `a2`. They also shape the filter response and are converted to Q15 fixed-point before filtering. |

The input samples (`sample_input`) and the expected reference output (`expected_output`) are generated for the chosen `INPUT_SIZE` and coefficient set. The output is compared against the reference using an average absolute error threshold (the first 100 samples are skipped to let the filter settle). If you change the size, the input data, or the coefficients, the reference data must be regenerated to match the new correct result.

---

## 5. Measuring Energy on the EVK

`main.c` is instrumented to profile energy on the E1x EVK following [Measuring application energy on the E1x EVK](https://docs.efficient.computer/programming/measure-energy). It drives the `AON4` GPIO to mark the code region under measurement:

- `AON4` is held **high** during setup, then driven **low** immediately before the measured loop (falling edge = start).
- The `biquad_filter` kernel runs `NUM_ITERATIONS` times inside the region.
- `AON4` is driven **high** immediately after the loop (rising edge = end).

The correctness check and all `printf` output stay outside the measured region. The AON GPIO calls are guarded by `#ifdef __EFFCC__`, so the app still builds and runs on a host without the EVK SDK.

To take a measurement:

1. Capture the EVK data serial port to a CSV file (using the udev rule device name if installed):

   ```bash
   minicom --device /dev/eff-power -C capture.csv
   ```

2. From another terminal, build and flash the app, then let it run to completion and close the serial client.

3. Process the capture with the same iteration count used by the firmware:

   ```bash
   python3 process-energy.py capture.csv --iters 1000
   ```

`process-energy.py` reports the runtime and energy per iteration on the E1x chip rail (`power_e1x(mW)`) by default; pass `--energy-column` to select a different rail (see the rail table in the measurement guide). If the reported runtime is not comfortably longer than the EVK sample period, raise `NUM_ITERATIONS` in `main.c` and pass the new value to `--iters`.
