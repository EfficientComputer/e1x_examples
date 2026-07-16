# Levenshtein Edit Distance (Levenshtein)

This example computes the **Levenshtein edit distance** between two strings on the Electron E1x general-purpose processor. It measures how many single-character edits separate the strings, showing how a table-filling dynamic program maps onto the Fabric architecture.

---

## 1. Overview

### What is edit distance?

The Levenshtein edit distance between two strings is the smallest number of single-character edits (insertions, deletions, or substitutions) needed to turn one string into the other. It is a standard measure of how similar two pieces of text are.

The computation fills a table where entry `dist[i][j]` holds the edit distance between the first `i` characters of one string and the first `j` characters of the other. Each entry is derived from its three neighbors (above, left, and above-left), so the full distance is built up from the distances of shorter prefixes.

### Mathematical Definition

    dist[i][j] = min(
        dist[i-1][j]   + 1,             (deletion)
        dist[i][j-1]   + 1,             (insertion)
        dist[i-1][j-1] + cost           (match or substitution)
    )

Where:

- `cost` is 0 if the two characters being compared are equal, and 1 otherwise
- the first row and first column are initialized to `0, 1, 2, ...`, the cost of building a prefix from an empty string
- the answer is the bottom-right entry, `dist[len1][len2]`

---

## 2. Why This Kernel Matters

Edit distance is a core building block in text and sequence processing:

- **Spelling correction and autocomplete**, where the closest dictionary word to a typo is chosen
- **Search and fuzzy matching**, where near-matches are ranked by similarity
- **Bioinformatics**, where DNA and protein sequences are aligned
- **Data cleaning and deduplication**, where near-duplicate records are merged
- **Version comparison**, where the differences between two texts are measured

It is a representative dynamic program whose work is dominated by comparisons and running minimums over a two-dimensional table rather than dense floating-point arithmetic.

---

## 3. Why EFF Hardware Performs Well

The Electron E1x runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

The edit-distance table has a regular dependency structure that maps well onto the Fabric:

- Entries along a diagonal depend only on earlier diagonals, so they can be computed **in parallel**
- The character comparisons and the minimum selections **overlap**, so control and memory work proceed together
- Recently computed table values stay **resident** near the compute elements, cutting repeated trips to memory
- The take-the-minimum-of-three logic at each entry flows as a **pipeline**, keeping control overhead low

The result is steady progress on a comparison-heavy dynamic program at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## 4. Configurable Parameters

These definitions in `main.c` control the benchmark. Change them to resize the problem or re-run it.

| Definition       | Default | Effect                                                                                                                                                                      |
| ---------------- | ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `NUM_ITERATIONS` | `1`     | How many times the kernel runs. Increase it to average out measurement noise when benchmarking.                                                                             |
| `len1`           | `28`    | The length of the first string. It must match the actual length of the `s1` literal.                                                                                        |
| `len2`           | `27`    | The length of the second string. It must match the actual length of the `s2` literal.                                                                                       |
| `EDIT_DISTANCE`  | `5`     | The expected edit distance the result is compared against for the pass/fail print. If you change the input strings, this must be updated to match the new correct distance. |

The two input strings are hardcoded literals in `main.c`. You can edit these literals to compare different text, but if you do, update `len1` and `len2` to the new string lengths and update `EDIT_DISTANCE` to the new correct result.
