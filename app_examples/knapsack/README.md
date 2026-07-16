# 0/1 Knapsack (Knapsack)

This example solves the **knapsack problem** with dynamic programming on the Electron E1x general-purpose processor. It selects items to maximize total value without exceeding a weight capacity, showing how a table-filling dynamic program maps onto the Fabric architecture.

---

## 1. Overview

### What is the knapsack problem?

The knapsack problem asks which subset of items to pack so that their combined value is as large as possible while their combined weight stays within a fixed capacity. Each item has a weight and a value. In this example each item also has a limited count of copies available, so a bounded number of each item may be taken.

The solution uses dynamic programming: it builds a table where each entry holds the best value achievable using the first `i` items within a given weight budget `j`. Every entry is computed from entries already filled in, so the final answer is assembled from smaller subproblems that are each solved only once.

### Mathematical Definition

    m[i][j] = max over k of ( m[i-1][j - k * weight[i]] + k * value[i] )

Where:

- `m[i][j]` is the best total value using the first `i` items within weight budget `j`
- `k` ranges from 0 up to the available count of item `i`, subject to `k * weight[i] <= j`
- `k = 0` reduces to `m[i-1][j]`, meaning item `i` is not taken

After the table is filled, the chosen counts are recovered by walking backward from the final entry and comparing each row to the row above it.

---

## 2. Why This Kernel Matters

Knapsack-style optimization appears throughout resource allocation and decision-making:

- **Budgeting and investment**, where a fixed budget is spread across projects to maximize return
- **Cargo and storage loading**, where limited capacity is packed for maximum value
- **Resource scheduling**, where tasks are chosen to fit within a time or memory limit
- **Cutting and packing**, where material is allocated to jobs to minimize waste
- **Feature and configuration selection**, where a value target is met under a cost constraint

It is a good example of a data-dependent dynamic program, where the work is dominated by table lookups, comparisons, and running maximums rather than dense floating-point arithmetic.

---

## 3. Why EFF Hardware Performs Well

The Electron E1x runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

The knapsack table is filled row by row, with many independent entries per row, and that structure maps well onto the Fabric:

- The entries across a single capacity row are computed **in parallel**, since each depends only on the previous row
- The running-maximum comparisons and the table lookups **overlap**, so memory access is hidden behind useful work
- The previous row of the table stays **resident** near the compute elements, reducing repeated trips to memory
- The compare-and-select logic that keeps the best value flows as a **pipeline**, keeping control overhead low

The result is steady progress on a comparison-heavy dynamic program at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## 4. Configurable Parameters

These definitions in `main.c` control the benchmark. Change them to resize the problem or re-run it.

| Definition       | Default | Effect                                                                                                                                           |
| ---------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| `NUM_ITERATIONS` | `1`     | How many times the solver runs. Increase it to average out measurement noise when benchmarking.                                                  |
| `n`              | `5`     | The number of items considered. This sets the number of rows in the dynamic-programming table.                                                   |
| `w`              | `400`   | The weight capacity of the knapsack. This sets the number of columns in the table, so a larger value means a larger table and more work per run. |

The items themselves (each item's weight, value, and available count) are hardcoded literals in the `init_items` function in `main.c`. You can edit these literals to change the problem. The correctness check compares the result against an expected total count of `6`, total weight of `395`, and total value of `730`, all hardcoded in `main.c`. If you change the items, the capacity, or the item count, these expected values must be updated to match the new correct answer.
