# Merge Sort (Sort)

This example sorts an array of 32-bit integers into ascending order using a bottom-up **merge sort** on the Electron E1x general-purpose processor. It repeatedly merges sorted runs of doubling length, showing how a comparison-heavy kernel dominated by data movement maps onto the Fabric architecture.

---

## Overview

### What is merge sort?

Merge sort builds a sorted array by repeatedly merging runs that are already sorted. Two sorted runs can be combined into one longer sorted run in a single pass: compare the front element of each run, take the smaller one, and advance only the side it came from. Starting from single elements, which are trivially sorted, the run length doubles on every pass until the whole array is one sorted run.

This example uses the bottom-up form of the algorithm, so the passes are driven by loops rather than by recursion. It merges between two buffers. Each pass reads runs from one buffer and writes merged runs into the other, then the two buffers swap roles for the next pass. Because output is never written over input the same pass is still reading, no element has to be shifted inside the array. If the last pass leaves the sorted data in the scratch buffer, it is copied back into the original array.

### Procedure

1. Treat each single element as a sorted run of length 1.
2. For the current run length `k`, walk the array in slices of `2k` elements, pairing each run with the run that follows it.
3. Merge each pair by comparing the two front elements, writing the smaller one to the output buffer, and advancing only the side that was taken.
4. Once one run is exhausted, copy the remaining elements of the other run straight across.
5. Swap the two buffers, double `k`, and repeat until `k` spans the whole array.

For `n` elements, every pass touches all `n` elements and the run length doubles, so there are `log2(n)` passes and the total work is on the order of `n log n` comparisons. That bound holds for every input, including the worst case.

---

## Why This Kernel Matters

Sorting is one of the most common operations in computing, and it underpins many other algorithms:

- **Databases and search**, where sorted order makes lookups, ranking, and merge joins fast
- **Sensor and event processing**, where readings are ordered by time or magnitude before analysis
- **Statistics**, where medians, percentiles, and other order statistics need ranked data
- **Graphics and geometry**, where primitives are ordered by depth or coordinate before rendering
- **Compression and data pipelines**, where grouping equal or nearby values first makes later stages more effective

It is a good measure of how well an architecture handles data movement and data-dependent control, because the work is comparisons, index updates, and copies rather than dense arithmetic.

---

## Why Efficient Hardware Performs Well

The Electron E1x runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

Merge sort is dominated by comparisons, index updates, and memory traffic rather than heavy arithmetic, and the Fabric handles that mix well:

- The slices merged within a single pass are independent of one another, so many of them advance **at the same time** across the Fabric.
- The inner merge loop is written without branches, so the comparison result selects both the value that is stored and the index that advances, and the merge flows as a **pipeline** instead of stalling on an unpredictable decision.
- Loading the next candidate element **overlaps** with writing the current output, so memory latency is hidden behind useful work.
- The two front values and the running indexes stay **resident** near the compute elements, cutting repeated trips to memory.
- Alternating between the two buffers keeps the access pattern of each pass sequential and predictable, so little energy goes to control overhead.

The result is efficient sorting at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## Configurable Parameters

These definitions in `main.c` control the benchmark. Change them to resize the problem or re-run it.

| Definition       | Default | Effect                                                                                                                                                                              |
| ---------------- | ------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `NUM_ITERATIONS` | `1`     | This is how many times the sort runs. Increase it to average out measurement noise when benchmarking.                                                                                |
| `INPUT_SIZE`     | `512`   | This is the number of elements sorted. This sets the problem size: a larger value means more elements per pass and more passes. It also sizes the input, working, and scratch arrays. |

