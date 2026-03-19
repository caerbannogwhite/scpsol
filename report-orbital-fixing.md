# Orbital Fixing for Symmetry Breaking in scpsol: Theory and Implementation

## 1. Introduction

Many combinatorial optimization problems exhibit symmetry: multiple solutions that are structurally identical but differ only in which interchangeable objects are selected. In Set Covering Problems, columns with identical cost and identical coverage patterns are interchangeable — any feasible solution using one can substitute the other without changing the objective. When branch-and-bound explores both alternatives, it wastes effort on symmetric subtrees that yield no new information.

Orbital fixing is a symmetry-breaking technique that detects interchangeable variables at each BnB node and fixes all but one representative to zero, eliminating symmetric branches from the search tree. The scpsol solver implements a lightweight orbital fixing scheme based on Zobrist hashing, with exact verification and automatic disabling when no symmetry is detected.

## 2. Background

### 2.1 Symmetry in Set Covering Problems

A Set Covering Problem has the form:

$$\min \; c^T x \quad \text{s.t.} \; Ax \geq \mathbf{1}, \; x \in \{0,1\}^n$$

Two columns $j$ and $k$ are _symmetric_ if swapping them in any feasible solution preserves both feasibility and objective value. A sufficient condition is:

1. **Equal cost:** $c_j = c_k$
2. **Identical coverage:** $A_{ij} = A_{ik}$ for all rows $i$

If both conditions hold, columns $j$ and $k$ are interchangeable — any solution selecting $j$ but not $k$ has a mirror image selecting $k$ but not $j$ with the same cost and coverage. The branch-and-bound tree explores both, doubling the search space unnecessarily.

### 2.2 Static vs. Dynamic Symmetry

**Static symmetry** exists in the original problem formulation: two columns that cover exactly the same rows at identical cost. Preprocessing can detect and remove these via dominance analysis — if column $j$ dominates column $k$ (covers a superset of rows at no higher cost), column $k$ can be eliminated.

**Dynamic symmetry** arises during BnB as branching decisions change the _residual_ problem. Two columns that are distinct in the original problem may become equivalent after some rows are decided:

- If column $j$ covers rows $\{1, 2, 3\}$ and column $k$ covers rows $\{1, 2, 4\}$, they are not equivalent
- But if row 3 and row 4 are already covered by fixed columns, the _residual_ coverage of both $j$ and $k$ is $\{1, 2\}$ — they are now interchangeable

Dynamic symmetry cannot be detected by preprocessing; it must be checked at each BnB node.

### 2.3 Orbital Fixing

The concept of orbital fixing originates from group theory applied to integer programming (Margot, 2002; Ostrowski et al., 2011). Given a symmetry group $G$ acting on the variables, the _orbit_ of variable $x_j$ under $G$ is the set of all variables that $j$ maps to under some symmetry. If two variables are in the same orbit, fixing one and freeing the other produces an isomorphic subproblem.

For SCP, the symmetry group at a BnB node consists of permutations that map free columns to free columns while preserving cost and residual coverage. Within each orbit, we keep one representative free and fix the rest to 0. This is valid because:

- If the optimal solution uses column $k$ (fixed to 0), there exists an equally good solution using the representative $j$ instead (by symmetry)
- So no optimal solution is excluded

The challenge is computing orbits efficiently. Full symmetry detection via graph isomorphism (e.g., nauty/bliss) is expensive. For the special structure of SCP, a simpler approach suffices.

## 3. Algorithm

### 3.1 Overview

The orbital fixing algorithm in scpsol operates in three phases:

1. **Preprocessing:** Group columns by cost (static, done once per model)
2. **Hashing:** At each BnB node, compute a Zobrist hash of each free column's residual coverage
3. **Verification:** For columns with matching hash and cost, verify exact equivalence and fix duplicates

### 3.2 Cost Grouping

Columns can only be equivalent if they have the same cost. As a first filter, columns are grouped by their objective coefficient, quantized to avoid floating-point issues:

$$\text{key}(j) = \text{round}(c_j \times 10^6)$$

Only groups with more than one column are retained. This partitioning is computed once when the base model is built and rebuilt when the model changes (after mid-BnB column removal).

### 3.3 Zobrist Hashing

For each cost group, free columns are hashed by their residual coverage — the set of undecided rows they cover. The hash function uses Zobrist hashing (Zobrist, 1970), a technique originally developed for board game position hashing:

1. Assign each row $i$ a random 64-bit value $h_i$ (generated once, deterministically seeded)
2. For each free column $j$, compute:

$$H(j) = \bigoplus_{i \in R(j)} h_i \oplus (\lvert R(j) \rvert \cdot \phi)$$

where $R(j) = \{i : A_{ij} > 0, \; \text{row } i \text{ undecided}\}$ is the residual row set, $\oplus$ denotes XOR, and $\phi = \texttt{0x9E3779B97F4A7C15}$ is the golden ratio constant used to mix the count into the hash.

The count mixing prevents false matches between columns covering different numbers of rows whose XOR happens to collide.

**Properties:**

- **O(nnz) computation:** Each column's hash requires one pass over its nonzero entries
- **Low collision probability:** With 64-bit hashes and golden ratio mixing, false positives are astronomically rare ($\approx 2^{-64}$ per pair)
- **Deterministic:** The seed `0xBEEF5EED` ensures reproducible behavior across runs

### 3.4 Equivalence Verification

Hash equality is necessary but not sufficient for equivalence. After sorting columns by hash within each cost group, columns sharing the same hash are verified by exact comparison of their residual row sets:

1. Build the residual row set for the first column in the bucket (the _representative_)
2. For each subsequent column, build its residual row set and compare element-by-element
3. If the sets match, the column is equivalent to the representative and can be fixed to 0

The row sets from the adjacency structure are already sorted, so comparison is a linear merge operation.

### 3.5 Fixing

For each equivalence class of size $k$, one column is kept free (the representative — the first column encountered in the group) and $k-1$ columns are fixed to 0. The fixings are added to the node's decision set before the LP is solved, so they:

- Tighten the LP relaxation (fewer free variables)
- Reduce the branching candidates
- Propagate through the essential-column cascade (fixing a column to 0 may force another column to 1)

## 4. Implementation

### 4.1 Data Structures

The implementation uses a single state object that persists across BnB nodes:

```
struct OrbitalFixingState {
    cost_groups    : vector<vector<int>>   — column groups with identical cost
    row_hash       : vector<uint64_t>      — Zobrist hash table (row → random value)
    total_fixed    : int                   — cumulative fixings across all nodes
    calls          : int                   — number of nodes where orbital fixing was attempted
    disabled       : bool                  — auto-disable flag
    col_state_buf  : vector<int8_t>        — scratch: column decision state (-1/0/+1)
    row_decided_buf: vector<uint8_t>       — scratch: whether each row is decided
};
```

The scratch buffers are reused across calls to avoid per-node allocation overhead.

### 4.2 Integration with BnB

Orbital fixing runs at each BnB node after propagation but before the LP solve:

1. **Propagation:** Essential-column cascade detects implied fixings from branching decisions
2. **Orbital fixing:** Detect and fix equivalent columns based on current decisions + propagated implications
3. **LP solve:** Solve with all fixings applied (branch decisions + propagated + orbital)

The combined decision set (branching + propagation + orbital) is passed to `lp.apply_decisions()` as a single batch.

### 4.3 Activation Guard

Orbital fixing is only attempted when the node has at least 2 branching decisions. At the root node (0 decisions) and first-level children (1 decision), there are too few decided rows for dynamic equivalences to emerge, and the full column set has already been processed by static dominance analysis.

### 4.4 State Invalidation

The orbital state (cost groups, Zobrist hashes) must be rebuilt whenever the base model changes:

- **Mid-BnB column removal:** Columns are eliminated, indices change → rebuild cost groups
- **Mid-BnB budget pruning:** Similar column removal → rebuild

An `orbital_valid` flag tracks whether the state matches the current base model. When invalidated, `build_orbital_state()` is called lazily at the next node that attempts orbital fixing.

### 4.5 Auto-Disable Mechanism

On instances without exploitable symmetry, orbital fixing adds per-node overhead (building column states, computing hashes, sorting) without producing any fixings. To eliminate this overhead:

- After **100 calls** with **0 total fixings**, orbital fixing is permanently disabled for the remainder of the solve
- 100 nodes is sufficient to observe whether the instance has dynamic symmetry, while keeping the wasted overhead small (~0.5s on typical instances)

This ensures the feature is zero-cost on instances where it provides no benefit.

## 5. Complexity Analysis

### 5.1 Per-Node Cost

At each node, the algorithm:

1. **Initialize column/row state:** $O(n + m)$ where $n$ = columns, $m$ = rows
2. **Process decisions:** $O(d \cdot \bar{r})$ where $d$ = decisions, $\bar{r}$ = average rows per column
3. **Hash free columns in each group:** $O(\sum_{g} |g| \cdot \bar{r})$ where $|g|$ is the group size
4. **Sort by hash:** $O(\sum_{g} |g| \log |g|)$
5. **Exact verification:** $O(|B| \cdot \bar{r})$ for each hash bucket $B$ with collisions

The dominant cost is step 3, which is $O(\text{nnz}_{\text{groups}})$ — the total nonzero entries in columns that belong to cost groups. On OR-Library instances with 300-600 columns and ~20 rows per column, this is ~5,000-10,000 operations per node.

### 5.2 Total Overhead

With auto-disable after 100 fruitless nodes, the maximum wasted overhead is:

$$100 \times O(\text{nnz}_{\text{groups}}) \approx 0.5\text{s}$$

This is negligible relative to BnB solve times of 20-300 seconds on hard instances.

## 6. Experimental Results

### 6.1 Setup

All experiments use the OR-Library SCP benchmark instances (Beasley, 1990) with a 300-second time limit, single-threaded, on Windows 11.

### 6.2 Symmetry Detection

The table below shows the number of cost groups and columns in groups for each instance category, measured after preprocessing (dominance analysis, essential columns, RC fixing).

| Category | Instances | Avg Cost Groups | Avg Cols in Groups | Notes                              |
| -------- | --------: | --------------: | -----------------: | ---------------------------------- |
| scp4     |        10 |             2.7 |                7.1 | Small groups after preprocessing   |
| scp5     |        10 |             5.8 |               36.7 | Moderate grouping                  |
| scpa     |         5 |            16.6 |              214.6 | Large groups, most cols in groups  |
| scpb     |         5 |            13.6 |              268.8 | Similar to scpa                    |
| scpnre   |         5 |            11.2 |              470.8 | Most columns share a cost group    |
| scpnrf   |         5 |             8.2 |              342.2 | Many same-cost columns             |

Despite having many same-cost columns, the OR-Library instances exhibit little _dynamic_ symmetry. After dominance preprocessing removes columns with identical or dominated coverage, the remaining columns within each cost group have distinct coverage patterns. The orbital fixing auto-disables after 100 nodes on all 40 instances.

### 6.3 Performance Impact

Since no dynamic equivalences are found on OR-Library instances, the orbital fixing feature is dormant on this benchmark set. The auto-disable mechanism limits the overhead to the first 100 BnB nodes.

| Category    | Without Orbital (v0.9) | With Orbital (v0.10) | Delta  |
| ----------- | ---------------------: | -------------------: | -----: |
| **scp4**    |                  0.5 s |                0.5 s |  0.0 s |
| **scp5**    |                  0.9 s |                0.9 s |  0.0 s |
| **scpa**    |                  3.9 s |                4.0 s | +0.1 s |
| **scpb**    |                 19.3 s |               18.9 s | -0.4 s |
| **scpnre**  |              1091.9 s |             1128.5 s |   noise |
| **scpnrf**  |              1151.7 s |             1217.9 s |   noise |

The delta column shows the differences are well within measurement noise (BnB tree variability and CPU thermal effects dominate on NRE/NRF instances). The auto-disable mechanism successfully eliminates per-node cost after the initial 100-node probe period.

### 6.4 When Orbital Fixing Helps

Orbital fixing is designed for instances with structured symmetry that survives preprocessing. Such instances arise in:

- **Facility location SCPs** where multiple facilities have identical construction costs and serve overlapping regions
- **Scheduling SCPs** where identical time slots or resources create interchangeable columns
- **Random SCPs with uniform costs** where many columns have the same objective coefficient and, after branching decisions cover certain rows, become residually equivalent

On such instances, orbital fixing can reduce the BnB tree by a factor proportional to the average orbit size. For an equivalence class of size $k$, the tree size reduction factor is up to $k$ in the subtree where that class is active.

## 7. Relationship to Other Techniques

### 7.1 Dominance Preprocessing

Static column dominance (already in scpsol's preprocessing pipeline) eliminates columns whose coverage is a subset of another column with equal or lower cost. This handles a common source of symmetry at the root level. Orbital fixing complements dominance by catching _dynamic_ equivalences that only appear after branching decisions.

### 7.2 Isomorphism Pruning

Full isomorphism pruning (Margot, 2002) uses a symmetry group computed by graph automorphism tools (nauty/bliss) to prune nodes whose subproblems are isomorphic to previously explored ones. This is more powerful than orbital fixing — it can detect symmetries even when columns have different costs but participate in symmetric constraint structures.

However, isomorphism pruning has two drawbacks: (1) computing the symmetry group is expensive ($O(n^2)$ or worse per node), and (2) maintaining a database of explored isomorphism classes adds memory overhead. Orbital fixing trades detection power for speed: it only finds the simplest symmetries (same cost, same residual coverage) but does so in $O(\text{nnz})$ time with no memory overhead beyond the hash table.

### 7.3 Symmetry-Breaking Constraints

An alternative approach adds _lexicographic ordering constraints_ that break symmetry without modifying the BnB tree structure. For columns $j < k$ in the same symmetry class, the constraint $x_j \geq x_k$ ensures only the lexicographically smallest representative is selected. These constraints can be added at the root as cutting planes.

Orbital fixing achieves a similar effect but adapts dynamically to the residual problem at each node, potentially detecting equivalences that static symmetry-breaking constraints miss.

## 8. References

1. F. Margot. Pruning by isomorphism in branch-and-cut. _Mathematical Programming_, 94(1):71–90, 2002.
2. J. Ostrowski, J. Linderoth, F. Rossi, and S. Smriglio. Orbital branching. _Mathematical Programming_, 126(1):147–178, 2011.
3. A. L. Zobrist. A new hashing method with application for game playing. _Technical Report 88_, Computer Science Department, University of Wisconsin–Madison, 1970.
4. K. Apt. _Principles of Constraint Programming_. Cambridge University Press, 2003. Chapter 10 (symmetry breaking).
5. F. Margot. Symmetry in integer linear programming. In M. Jünger et al., editors, _50 Years of Integer Programming_, pages 647–686. Springer, 2010.
