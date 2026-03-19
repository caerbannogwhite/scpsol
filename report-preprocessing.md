# Preprocessing Techniques in scpsol: Theory and Implementation

## 1. Introduction

Preprocessing (or presolving) reduces the size of a combinatorial optimization problem before the main solver begins. For Set Covering Problems (SCP), this means removing redundant rows (constraints), dominated columns (variables), and fixing variables whose optimal values can be deduced. Effective preprocessing can dramatically reduce the search space and solve time, sometimes resolving the problem entirely without branch-and-bound.

The scpsol solver implements a multi-stage preprocessing pipeline that iterates multiple techniques until a fixpoint is reached — that is, until no further reductions are possible. This report describes the theory behind each technique, the implementation in scpsol, and how they interact.

## 2. Problem Formulation

The Set Covering Problem (SCP) is:

$$\min \; c^T x \quad \text{s.t.} \; Ax \geq \mathbf{1}, \; x \in \{0,1\}^n$$

where $A \in \{0,1\}^{m \times n}$ is the constraint matrix, $c > 0$ is the cost vector, and $\mathbf{1}$ is the all-ones vector. Each row $i$ must be covered by at least one selected column. We denote the set of columns covering row $i$ as $S_i = \{j : A_{ij} = 1\}$ and the set of rows covered by column $j$ as $R_j = \{i : A_{ij} = 1\}$.

## 3. Preprocessing Techniques

### 3.1 Essential Column Fixing

**Theory.** A column $j$ is _essential_ for row $i$ if it is the only column that can cover $i$, i.e., $S_i = \{j\}$. In any feasible solution, $x_j = 1$ must hold. Fixing $x_j = 1$ allows removing all rows in $R_j$ from the problem, since they are now guaranteed to be covered.

Crucially, removing rows can make other columns essential. If column $k$ was one of two columns covering row $i'$, and the other column was just fixed, then $k$ becomes essential for $i'$. This creates a _cascade_ of fixings.

**Implementation.** The essential column cascade in `row_reduce()` (`src/preprocessor.cpp`) iterates until quiescence. Each round scans all active rows, counts active covering columns, and fixes any column that is the sole cover for some row. Fixing a column removes all its rows, which may create new singleton-covered rows in the next round. The loop terminates when a full pass produces no new fixings. The fixed columns' costs are accumulated in `fixed_preprocess_cost` and subtracted from the incumbent bound, since they are known to be in every optimal solution.

### 3.2 Row Domination

**Theory.** Row $i$ _dominates_ row $i'$ if every column covering $i$ also covers $i'$, i.e., $S_i \subseteq S_{i'}$. Any feasible solution that satisfies the constraint for $i$ automatically satisfies the constraint for $i'$, so $i'$ is redundant and can be removed.

**Implementation.** In `row_reduce()`, active rows are sorted by coverage size (ascending), and for each pair $(i, i')$ where $|S_i| \leq |S_{i'}|$, a bitset subset check determines if $S_i \subseteq S_{i'}$. If so, $i'$ is deactivated. The bitset representation uses $\lceil n/64 \rceil$ words per row (where $n$ is the number of columns), enabling fast subset checks even for large instances.

### 3.3 Probing (Constraint Propagation)

**Theory.** For each active column $j$, we temporarily assume $x_j = 0$ and propagate the consequences via the essential column cascade. If this leads to a row with no remaining covering columns (infeasibility), then $x_j = 0$ is impossible in any feasible solution, so $x_j = 1$ must hold. This is equivalent to a single round of _node presolve_ or _probing_ in MIP solvers.

More formally: setting $x_j = 0$ may force other columns to be essential (creating $x_k = 1$), which removes rows and may make further columns essential. If the cascade reaches a contradiction (an uncoverable row), the original assumption $x_j = 0$ is refuted.

**Implementation.** In `row_reduce()`, for each active column $j$:

1. Temporarily deactivate $j$
2. Run the essential column cascade, tracking all changes on a rollback stack
3. If any row becomes uncoverable, mark $j$ as fixed to 1
4. Rollback all temporary changes

After all probing, a final essential column cascade picks up any newly created essentials.

### 3.4 Single Column Dominance

**Theory.** Column $j$ is _dominated_ by column $k$ if $c_k \leq c_j$ and $R_j \subseteq R_k$ — that is, $k$ covers all the rows that $j$ covers, at no greater cost. In any optimal solution that uses $j$, replacing $j$ with $k$ yields a solution that is at least as good. Therefore $j$ can be removed.

**Implementation.** The `SingleColumnDominanceRule` checks every pair of active columns. For each target column $j$, if any candidate $k$ has $c_k \leq c_j$ and $R_j \subseteq R_k$, then $j$ is deactivated. Tie-breaking: if costs are equal, the lower-indexed column is kept.

**Bitset acceleration.** Subset checks are accelerated using `DynBitset` (`src/bitset_util.h`): each column's row-coverage is encoded as a fixed-width bitset of $\lceil m/64 \rceil$ words. The subset check $R_j \subseteq R_k$ becomes $\lceil m/64 \rceil$ bitwise AND+CMP operations, replacing the $O(|R_j|)$ sorted-merge scan. For $m = 300$ (scpb instances), this is just 5 word comparisons.

### 3.5 Two-Column Dominance

**Theory.** Column $j$ is _pair-dominated_ if there exist columns $k_1, k_2$ such that $c_{k_1} + c_{k_2} < c_j$ and $R_j \subseteq R_{k_1} \cup R_{k_2}$. Any solution using $j$ can be improved by using $k_1$ and $k_2$ instead, so $j$ is redundant.

**Implementation.** The `TwoColumnDominanceRule` iterates over all target columns and all pairs of candidates. Time limit checks prevent excessive computation on large instances.

**Bitset acceleration.** The union-covers check $R_j \subseteq R_{k_1} \cup R_{k_2}$ is performed using `DynBitset::is_subset_of_union()`, which verifies the condition using $\lceil m/64 \rceil$ bitwise OR, AND, and compare operations per pair.

### 3.6 Cost-Driven Replacement

**Theory.** This generalizes pair dominance by considering 2-column and 3-column replacements, focused on expensive columns first. By processing columns in descending cost order, the most impactful reductions are found first. A column $j$ is removed if any combination of 2 or 3 cheaper columns covers all its rows at lower total cost.

**Implementation.** The `CostDrivenReplacementRule` builds a reverse index (columns by row) to efficiently find candidate replacements. For each target column (processed from most expensive to cheapest), it collects all columns sharing at least one row with the target, sorts them by cost, and checks pairs and triples.

**Bitset acceleration.** Pair and triple union-covers checks use `DynBitset::is_subset_of_union()` and `DynBitset::is_subset_of_union3()`, which compute the union-subset test in $\lceil m/64 \rceil$ bitwise operations regardless of column cardinalities.

### 3.7 Greedy Multi-Column Dominance Finder

**Theory.** The pairwise (§3.4–3.6) dominance checks are exact but have $O(n^2)$ or $O(n^3)$ complexity, which becomes prohibitive for large column sets. The greedy multi-column dominance finder generalizes these checks to arbitrary numbers of replacement columns using a fast heuristic: for each candidate column $j$, it tests whether a set of cheaper "kept" columns $K$ can collectively cover all rows of $j$ at total cost $\leq c_j$.

The algorithm follows the structure of Algorithm 1 from Grossman and Wool (1997):

1. Initialize $K$ with all unit-cost columns ($c_j = 1$).
2. Sort remaining columns $N$ by ascending cost.
3. For each $j \in N$ (cheapest first):
   - Test: can columns in $K$ cover all rows $R_j$ at cost $\leq c_j$?
   - If yes: remove $j$ (it is dominated).
   - If no: add $j$ to $K$ (it is needed for future tests).

**Implementation.** The `dominance_finder()` function in `src/preprocessor.cpp` uses a greedy set cover heuristic for each dominance test, rather than solving a sub-MIP. For each candidate $j$:

1. **Feasibility check**: Verify every row of $j$ has at least one covering column in $K$. If not, $j$ cannot be dominated — keep it.
2. **Greedy cover**: Collect all $K$-columns sharing rows with $j$. Repeatedly select the column with the best (uncovered rows of $j$) / cost ratio until all rows are covered or the budget is exhausted.
3. If all rows are covered within budget $c_j$, remove $j$.

**Key implementation details:**

- Maintains a reverse index `tilde_cols_by_row` (updated incrementally as $K$ grows).
- Uses reusable `uncovered` and `is_candidate` buffers to avoid per-column allocation.
- Worst-case $O(|N| \cdot |K| \cdot \bar{r})$ where $\bar{r}$ is the average rows per column, but very fast in practice (under 0.1s on 3000-column instances).

**Impact.** On OR-Library benchmark instances, the greedy DF removes 50–85% of surviving columns (after cost/budget pruning) with negligible overhead. By running before the $O(n^2)$ pairwise checks, it reduces their input from $\sim$2600 to $\sim$400 columns, yielding a $\sim$5× speedup in total preprocessing time.

### 3.8 Incumbent Budget Pruning

**Theory.** Given a known upper bound $z^*$ (incumbent), any column $j$ with $c_j \geq z^*$ cannot appear in any improving solution, since it alone would exceed the budget. More generally, if including $j$ requires a minimum additional cost to cover the remaining rows that exceeds $z^* - 1$, then $j$ can be pruned.

The budget for column $j$ is $B_j = \lfloor z^* \rfloor - 1 - \lfloor c_j \rfloor$. If $B_j < 0$, the column is immediately pruned. For small budgets ($B_j < 1$), we check whether a single cost-1 column can cover all rows not covered by $j$. For larger budgets, we compute a lower bound on the cost of covering the remaining rows and prune if it exceeds $B_j$.

**Implementation.** The `IncumbentBudgetPruningRule` processes columns from most expensive to cheapest, building a per-row sorted column list for efficient minimum-cost lookups.

### 3.9 Reduced-Cost Fixing

**Theory.** After solving the LP relaxation with dual bound $z_{LP}$ and incumbent $z^*$, the _gap_ is $\Delta = z^* - z_{LP}$. For any column $j$ with LP value $x_j = 0$ and reduced cost $\bar{c}_j > \Delta$, fixing $x_j = 1$ would increase the LP bound beyond $z^*$, so $x_j = 0$ in every optimal solution. This is a standard technique in mixed-integer programming.

**Implementation.** In `solver.cpp` (Phase 4-5), reduced costs are computed from the LP dual solution. Columns at their lower bound ($x_j \approx 0$) with $\bar{c}_j > \Delta + \epsilon$ are removed.

### 3.10 Cost Reduction

**Theory.** Any column with $c_j \geq z^*$ (the incumbent bound) can be removed, since including it alone would match or exceed the best known solution. Combined with the SCP constraint that at least one other column must be present, the total would strictly exceed $z^*$.

**Implementation.** The `apply_cost_reduction()` function removes all columns with cost at or above the incumbent bound.

## 4. Greedy Heuristic

### 4.1 Theory

The greedy algorithm by Chvátal (1979) provides an initial feasible solution and upper bound. At each iteration, it selects the column with the best _cost-effectiveness ratio_:

$$j^* = \arg\max_{j} \frac{|\{i \in R_j : i \text{ uncovered}\}|}{c_j}$$

That is, it picks the column that covers the most uncovered rows per unit cost. This has a well-known approximation guarantee: the greedy solution costs at most $H(m)$ times the optimal, where $H(m) = 1 + 1/2 + \cdots + 1/m$ is the $m$-th harmonic number.

### 4.2 Implementation

The `greedy_set_cover_heuristic()` function in `src/preprocessor.cpp` implements the Chvátal greedy with **incremental coverage tracking**:

1. Build a reverse index (rows covered by each column)
2. Initialize per-column coverage counts `cov_count[j]` = number of uncovered rows column $j$ covers
3. At each step, scan all unused columns, select the one with the highest `cov_count[j] / c_j` ratio
4. Mark covered rows, and for each newly covered row, decrement `cov_count[k]` for all columns $k$ covering that row
5. Repeat until all rows are covered

The incremental update (step 4) ensures that coverage counts are maintained in $O(\text{nnz})$ total across all iterations, rather than recomputing from scratch at each step.

**Post-processing: Redundancy Removal.** After the greedy completes, a redundancy removal pass tries to remove unnecessary columns. Selected columns are sorted by cost (most expensive first), and each column is tentatively removed: if all rows remain covered (each row's coverage count stays $\geq 1$), the column is permanently dropped. This reduces the initial incumbent $z^*$, which in turn enables more aggressive cost reduction and budget pruning in subsequent preprocessing stages.

### 4.3 BnB Heuristic: Dual-Guided Cover Repair

**Theory.**

The `DualGuidedCoverRepairHeuristic` (`src/heuristics.cpp`) constructs integer-feasible solutions during branch-and-bound by combining LP relaxation information with greedy repair. Starting from the LP solution, it fixes near-integer variables and then greedily covers remaining rows using a score that combines uncovered-row gain and dual value gain, weighted by column cost.

**Implementation.** The heuristic uses a **column-to-rows transpose** (built once per invocation) and **incremental coverage tracking** for both the repair phase and the post-processing redundancy removal. When a column is added to the solution, only its covered rows' coverage values are updated; when a column is tentatively removed in post-processing, only its entries are subtracted. This reduces per-invocation cost from $O(n \times m \times w)$ (where $w$ is average row width) to $O(\text{nnz})$.

## 5. Preprocessing Pipeline

### 5.1 Iterated Fixpoint

A key design principle in scpsol is that preprocessing techniques are _iterated until fixpoint_. Each technique can enable further reductions by others:

- **Essential columns → dominance**: Fixing essential columns and removing their rows can make previously non-dominated columns now dominated.
- **Dominance → essential columns**: Removing dominated columns can leave rows with only one covering column, creating new essentials.
- **Probing → dominance**: Probing may fix columns and remove rows, enabling new dominance relationships.
- **Budget pruning → essential columns**: Removing expensive columns can create new essential columns.

The preprocessing loop runs all techniques in sequence and repeats the entire sequence until no further reductions occur (or a maximum of 10 rounds).

### 5.2 Pipeline Stages

The full preprocessing pipeline in `solver.cpp` proceeds as follows:

**Phase 1: Greedy Heuristic**

- Compute initial feasible solution and upper bound $z^*$
- Post-process with redundancy removal (remove unnecessary columns, most expensive first)

**Phase 2: Pre-LP Reduction (iterated)**
Repeat until fixpoint:

1. Cost reduction (remove columns with $c_j \geq z^*$)
2. Incumbent budget pruning
3. **Greedy multi-column dominance finder** (see §3.7)
4. Cost-driven replacement (2- and 3-column dominance)
5. Configured dominance rules (single, two-column)
6. Row reduction (essential columns, row domination, probing)

_Note:_ The greedy dominance finder (step 3) runs before the expensive pairwise/triplet checks (steps 4-5). This is critical for performance: it reduces the column count by 50-85% with negligible overhead, so the subsequent $O(n^2)$ and $O(n^3)$ rules operate on a much smaller problem.

**Phase 3: Root LP Relaxation**

- Solve LP relaxation on the reduced problem
- Try heuristics on root LP solution
- Record dual bound

**Phase 4-5: Post-LP Reduction (iterated)**
Repeat until fixpoint:

1. Reduced-cost fixing (using LP duals and gap)
2. Cost reduction
3. Incumbent budget pruning
4. Configured dominance rules
5. Row reduction

**Phase 6: Branch-and-Bound**

- Build final base model from reduced problem
- During BnB, mid-solve reductions triggered on incumbent improvement:
  - Cost-based column removal
  - Budget pruning with remapping of branch nodes

### 5.3 Mid-BnB Preprocessing

When a new incumbent is found during branch-and-bound, the tighter bound enables additional preprocessing. The functions `mid_bnb_column_removal()` and `mid_bnb_budget_pruning()` remove columns from the base model, remap all branch node decisions and cuts to the new variable indices, and prune frontier nodes whose dual bounds exceed the new incumbent.

## 6. Bitset Infrastructure

The `DynBitset` utility (`src/bitset_util.h`) provides a lightweight, variable-width bitset supporting subset and union-subset operations:

| Operation                      | Complexity              | Use                               |
| ------------------------------ | ----------------------- | --------------------------------- |
| `is_subset_of(other)`          | $O(\lceil n/64 \rceil)$ | Column dominance                  |
| `is_subset_of_union(a, b)`     | $O(\lceil n/64 \rceil)$ | Two-column dominance              |
| `is_subset_of_union3(a, b, c)` | $O(\lceil n/64 \rceil)$ | Cost-driven replacement (triples) |
| `popcount()`                   | $O(\lceil n/64 \rceil)$ | Coverage counting                 |

For standard benchmark instances with $m = 200$-$300$ rows, this means 4-5 uint64 word operations per check, replacing $O(|R_j|)$ sorted-merge operations. The `build_column_bitsets()` helper builds bitsets for all active columns from the `rows_by_column` transpose.

## 7. Node-Level Preprocessing

During branch-and-bound, scpsol applies lightweight preprocessing at each node:

### 7.1 Node-Level Propagation

Before solving a node's LP relaxation, the `propagate_fixings()` function (from `src/reliability.cpp`) runs an essential-column cascade on the node's fixed variables. This can:

- Detect infeasible nodes without an LP solve (when a row has no free columns)
- Discover implied fixings that tighten the LP

Propagation is only applied at nodes with $\geq 3$ decisions to avoid overhead on shallow nodes.

### 7.2 Node-Level Reduced Cost Fixing

After solving a node's LP, the solver computes reduced costs and fixes variables whose reduced cost exceeds the current gap ($z^* - z_{LP}^{\text{node}}$). These fixings are inherited by child nodes, progressively tightening the LP in deeper subtrees.

## 8. Configuration

The preprocessing behavior is controlled by the `SolverConfig` fields:

| Parameter               | Default        | Description                                  |
| ----------------------- | -------------- | -------------------------------------------- |
| `preprocess_rules`      | `"single,two"` | Comma-separated list of dominance rules      |
| `preprocess_time_limit` | `10.0`         | Time limit (seconds) per preprocessing phase |

Available rule tokens: `single` (single column dominance), `two` (two-column dominance), `cost_driven` (cost-driven replacement), `incumbent_budget` (budget pruning), `none` (disable all).

## 9. References

1. V. Chvátal. A greedy heuristic for the set-covering problem. _Mathematics of Operations Research_, 4(3):233–235, 1979.
2. E. Balas and A. Ho. Set covering algorithms using cutting planes, heuristics, and subgradient optimization. _Mathematical Programming Study_, 12:37–60, 1980.
3. T. Grossman and A. Wool. Computational experience with approximation algorithms for the set covering problem. _European Journal of Operational Research_, 101(1):81–92, 1997.
4. T. Achterberg, T. Koch, and A. Martin. Constraint integer programming: a new approach to integrate CP and MIP. _LNCS_, 3011:6–20, 2004.
5. M. Savelsbergh. Preprocessing and probing techniques for mixed integer programming problems. _ORSA Journal on Computing_, 6(4):445–454, 1994.
