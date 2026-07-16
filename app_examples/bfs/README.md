# Breadth-First Search (BFS)

This example performs a **breadth-first search** over a graph on the Electron E1x general-purpose processor. It labels every vertex with its discovery level (the number of hops from a chosen start vertex), showing how graph traversal maps onto the Fabric architecture.

---

## 1. Overview

### What is breadth-first search?

Breadth-first search explores a graph one level at a time. Starting from a single source vertex, it first visits all of the source's direct neighbors, then all of their unvisited neighbors, and so on outward. Because the search expands in rings of equal distance, the level at which a vertex is first reached equals its shortest hop distance from the source.

The graph is stored in Compressed Sparse Row (CSR) form using two arrays:

- an offsets array with one entry per vertex plus one, where entry `i` marks where vertex `i`'s neighbor list begins
- a neighbors array that concatenates every vertex's outgoing edge destinations end to end

To read the neighbors of vertex `i`, the traversal scans the neighbors array from `offsets[i]` up to `offsets[i + 1]`.

### Procedure

1. Set the start vertex's level to 0 and mark every other vertex as undiscovered.
2. For the current level, scan the vertices already discovered at that level.
3. For each such vertex, look at its neighbors. Any neighbor that is still undiscovered is assigned the next level.
4. Advance to the next level and repeat until no new vertices are discovered.

The output is the level array (the `frontier`), where each entry holds the discovery level of the corresponding vertex.

---

## 2. Why This Kernel Matters

Breadth-first search is a foundational graph operation that appears across many domains:

- **Shortest-hop routing**, where the fewest number of links between two points must be found
- **Network and social graph analysis**, where reachability and neighborhood structure are measured
- **Connectivity and clustering checks**, where components of a graph are identified
- **Planning and search**, where the closest reachable states are expanded first
- **Mesh and circuit analysis**, where influence spreads outward from a source node

It is a good stress test for irregular, data-dependent workloads because the work at each step depends on the graph structure rather than a fixed pattern.

---

## 3. Why EFF Hardware Performs Well

The Electron E1x runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

This suits breadth-first search, which is dominated by control decisions and memory access rather than heavy arithmetic:

- The scan of the frontier is split across the Fabric so that many vertices are examined **at the same time**
- Neighbor lookups and level updates **overlap**, so memory latency is hidden behind useful work
- The level array stays **resident** near the compute elements, reducing repeated trips to memory
- Testing whether a neighbor is undiscovered and assigning its level flows as a **pipeline**, keeping control overhead low

The result is steady progress on an irregular, branch-heavy traversal at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## 4. Configurable Parameters

These definitions control the benchmark. Change them to resize the problem or re-run it.

| Definition       | Default | Effect                                                                                                                                                                                 |
| ---------------- | ------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `NUM_ITERATIONS` | `1`     | How many times the traversal runs (in `main.c`). Increase it to average out measurement noise when benchmarking.                                                                       |
| `BFS_NUM_VTX`    | `1000`  | The number of vertices in the graph, defined in `bfs.h`. This sets the problem size. It is also used as the "undiscovered" marker in the level array, so the graph data must match it. |
| `BFS_NUM_EDGE`   | `10000` | The number of edges in the graph, defined in `bfs.h`. It sizes the neighbors array and must match the supplied graph data.                                                             |

The start vertex is a hardcoded literal in `main.c` (the traversal begins from vertex `2`). The graph itself (the `oa` offsets array and `na` neighbors array) is supplied from `bfs_graph.c`. The `final_frontier` array in `main.c` holds the expected discovery levels used for the pass/fail check. If you change the graph data, the start vertex, or the sizes, this expected result must be updated to match the new correct output.
