# Dijkstra's Shortest Path (Dijkstra)

This example computes the **shortest path** between two points on a weighted map using Dijkstra's algorithm on the Electron E1x general-purpose processor. It uses a heap-based priority queue and reports the path and its total distance, showing how a classic graph search maps onto the Fabric architecture.

---

## Overview

### What is Dijkstra's shortest path?

Dijkstra's algorithm finds the lowest-cost route from a start vertex to a destination vertex in a graph whose edges carry non-negative weights. It grows a set of vertices whose shortest distance from the start is already known, always expanding the closest not-yet-finalized vertex next. Because it never revisits a vertex once its distance is settled, the first time the destination is finalized its recorded cost is guaranteed to be the shortest.

In this example the graph is built from an ASCII map: characters mark walkable cells, `S` marks the start, `D` marks the destination, and each straight run between junctions becomes a weighted edge whose cost is the number of steps in the run.

### Procedure

1. Assign the start vertex a tentative distance of 0 and every other vertex an effectively infinite distance.
2. Keep candidate vertices in a priority queue (a min-heap) ordered by tentative distance.
3. Remove the vertex with the smallest tentative distance and mark its distance as final.
4. For each of its neighbors, offer a new candidate distance equal to the current vertex's distance plus the connecting edge's cost.
5. Repeat until the destination is finalized, then walk the recorded predecessors backward to recover the path.

A min-heap keeps step 3 efficient: inserting a candidate and removing the current minimum each cost time proportional to the logarithm of the number of entries, rather than a full scan.

---

## Why This Kernel Matters

Shortest-path search is one of the most widely used graph operations in real systems:

- **Navigation and routing**, where the lowest-cost route across a road or transit network is needed
- **Robot and drone path planning**, where an agent must reach a goal across a cost map while avoiding obstacles
- **Network routing**, where packets follow least-cost paths between nodes
- **Games and simulation**, where units move toward targets over varied terrain
- **Logistics and scheduling**, where least-cost sequences of steps are chosen

It is a good test of irregular, data-dependent work because the amount and order of computation depend on the graph and the edge weights rather than a fixed pattern.

---

## Why Efficient Hardware Performs Well

The Electron E1x runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

Dijkstra's algorithm is dominated by control decisions, comparisons, and pointer-driven memory access rather than heavy arithmetic, and the Fabric handles that mix well:

- The heap operations and the neighbor relaxations **overlap**, so comparison and memory work proceed together instead of stalling one behind the other
- Distances and predecessor links stay **resident** near the compute elements, cutting repeated trips to memory
- The compare-and-update logic that decides whether a neighbor's distance improves flows as a **pipeline**, keeping control overhead low
- Traversing each vertex's neighbor list maps directly onto the dataflow graph, so little energy goes to instruction handling

The result is steady progress on a branch-heavy, pointer-chasing search at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## Configurable Parameters

These definitions in `main.c` control the benchmark. Change them to resize the problem or re-run it.

| Definition       | Default | Effect                                                                                                                |
| ---------------- | ------- | --------------------------------------------------------------------------------------------------------------------- |
| `NUM_ITERATIONS` | `1`     | This is how many times the shortest-path search runs. Increase it to average out measurement noise when benchmarking.         |
| `MAP_WIDTH`      | `20`    | This is the width in characters of the ASCII map that defines the graph. It must match the actual layout of the `map` string. |
| `MAP_HEIGHT`     | `9`     | This is the height in rows of the ASCII map. It must match the actual layout of the `map` string.                             |

The map itself is a hardcoded string literal in `main.c`, with `S` marking the start and `D` marking the destination. You can edit this map to change the graph, but if you do, update `MAP_WIDTH` and `MAP_HEIGHT` to match. The correctness check compares the result against an expected path segment count of `8` and a total distance of `28`, both hardcoded in `main.c`. If you change the map, these expected values must be updated to match the new correct route.
