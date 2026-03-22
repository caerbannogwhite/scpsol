# Cut Generation in scpsol: Theory and Implementation

## 1. Introduction

Cutting planes (cuts) tighten the LP relaxation of a mixed-integer program by adding valid inequalities that are violated by the current fractional solution but satisfied by all integer-feasible solutions. For Set Covering Problems (SCP), cuts can significantly reduce the integrality gap at the root node, shrinking the branch-and-bound tree.

The scpsol solver implements three families of cutting planes: Chvátal-Gomory (CG) cuts derived from dual aggregation, row-pair CG cuts from constraint combinations, and Balas cover cuts from an additive branch-and-cut procedure. This report describes the theory behind each technique, the separation algorithms, and their integration into the solver.

## 2. Background

The Set Covering Problem (SCP) is $\min \; c^T x$ subject to $Ax \geq \mathbf{1}$, $x \in \{0,1\}^n$, where $A \in \{0,1\}^{m \times n}$ is the constraint matrix and $c > 0$ is the cost vector.

### 2.1 The LP Relaxation Gap

The SCP LP relaxation replaces $x_j \in \{0,1\}$ with $0 \leq x_j \leq 1$. The optimal LP value $z_{LP}$ provides a lower bound on the integer optimum $z^*$. The _integrality gap_ $z^* - z_{LP}$ determines the search effort: a tighter LP bound means fewer branch-and-bound nodes.

### 2.2 Valid Inequalities

A linear inequality $\alpha^T x \geq \beta$ is _valid_ for the SCP if every feasible integer solution satisfies it. A valid inequality is _violated_ by the current LP solution $\bar{x}$ if $\alpha^T \bar{x} < \beta$. Adding violated valid inequalities to the LP raises $z_{LP}$ without excluding any integer solution.

### 2.3 Reduced Costs

Given the LP dual solution $u \geq 0$, the reduced cost of variable $j$ is:

$$\bar{c}_j = c_j - \sum_{i : A_{ij} > 0} u_i A_{ij}$$

Variables with large positive reduced costs are expensive to include, while variables with reduced costs near zero are at or near their LP-optimal values. Reduced costs play a central role in the Balas additive procedure.

## 3. Chvátal-Gomory Cuts

### 3.1 Theory

The Chvátal-Gomory (CG) procedure derives valid inequalities by aggregating constraints with non-negative multipliers and rounding. Given the system $Ax \geq b$ with $x \geq 0$, and multipliers $u \geq 0$:

1. **Aggregate:** $\sum_i u_i (A^i x) \geq \sum_i u_i b_i$, yielding $\bar{a}^T x \geq \bar{b}$ where $\bar{a}_j = \sum_i u_i A_{ij}$
2. **Round coefficients:** Since $x_j \in \{0,1\}$, we can round each coefficient up: $\lceil \bar{a}_j \rceil x_j \geq \bar{a}_j x_j$
3. **Round RHS:** The resulting inequality $\sum_j \lceil \bar{a}_j \rceil x_j \geq \lceil \bar{b} \rceil$ is valid for all integer solutions

The cut is violated when $\sum_j \lceil \bar{a}_j \rceil \bar{x}_j < \lceil \bar{b} \rceil$. This can happen because the rounding of the RHS "jumps past" the current fractional value.

### 3.2 Dual Aggregated CG Separator

The `DualAggregatedCgSeparator` uses the LP dual solution as multipliers, which is a natural choice since dual variables capture the marginal value of each constraint.

**Algorithm:**

1. Select rows $i$ with dual value $u_i > \epsilon$
2. Compute aggregated coefficients: $\bar{a}_j = \sum_{i : u_i > \epsilon} u_i A_{ij}$ for all $j$
3. Compute aggregated RHS: $\bar{b} = \sum_{i : u_i > \epsilon} u_i b_i$
4. Check fractional part: $f_0 = \bar{b} - \lfloor \bar{b} \rfloor$. If $f_0 \notin (\epsilon, 1-\epsilon)$, skip (rounding won't help)
5. Round: coefficients $\lceil \bar{a}_j \rceil$ for nonzero terms, RHS $\lceil \bar{b} - \epsilon \rceil$
6. Check violation: $\sum_j \lceil \bar{a}_j \rceil \bar{x}_j < \lceil \bar{b} \rceil - \epsilon$
7. If violated, return the cut

This produces at most one cut per invocation, targeting the strongest single aggregation.

### 3.3 Row-Pair CG Separator

The `RowPairCgSeparator` considers pairwise combinations of constraints, exploring a broader space of aggregations than the single dual-weighted version.

**Algorithm:**

1. Identify rows with $u_i > \epsilon$ and sort by dual value (descending)
2. Limit to the top 40 rows to control computational cost
3. For each pair $(i_1, i_2)$ among the top rows:
   - Compute $\bar{b} = u_{i_1} b_{i_1} + u_{i_2} b_{i_2}$
   - Check fractional part condition
   - Aggregate coefficients from both rows
   - Round and check violation
4. Collect candidates, sort by violation magnitude (most violated first)
5. Return up to 30 cuts

The pairwise approach can find cuts invisible to the single aggregation because different constraint combinations expose different fractional structures.

## 4. Balas Cover Cuts

### 4.1 Theory

For the SCP, a _cover_ is a minimal set of columns $C$ such that every feasible solution must include at least one column from $C$. The corresponding _cover inequality_ is:

$$\sum_{j \in C} x_j \geq 1$$

This is valid because if all columns in $C$ were set to zero, at least one row would be uncovered. Cover inequalities are the SCP analogue of knapsack cover cuts in general MIP.

The challenge is finding covers that are _violated_ by the current LP solution — that is, sets $C$ where $\sum_{j \in C} \bar{x}_j < 1$.

### 4.2 The Balas Additive Procedure (BCG)

The Balas Cut Generation (BCG) algorithm, implemented in `balas_cut_generate()`, constructs violated cover cuts using an iterative additive procedure based on reduced costs. The algorithm is rooted in the work of Balas and Ho (1980) on set covering via cutting planes and subgradient optimization.

**Inputs:**

- LP primal and dual solutions
- Reduced costs $\bar{c}_j$
- Set $S$ of candidate variables: fractional ($\bar{x}_j > \epsilon$) with positive reduced cost ($\bar{c}_j > \epsilon$)
- Tight rows $T$ = binding constraints where $\sum_j A_{ij} \bar{x}_j \approx 1$
- Incumbent bound $z^*$

Let $N_i = \{j : A_{ij} = 1\}$ denote the set of columns covering row $i$.

**Algorithm:**

1. Initialize: working set $S$, cut set $W = \emptyset$, dual progress $y = 0$
2. While $S \neq \emptyset$ and $y < z^*$:
   1. Compute gap = $z^* - y$
   2. Find threshold $v_t$: the smallest reduced cost in $S$ that reaches the gap, i.e., $v_t = \min\{\bar{c}_j : j \in S,\; \bar{c}_j \geq \text{gap}\}$. If no column reaches the gap, use $v_t = \max\{\bar{c}_j : j \in S\}$
   3. Let $J = \{j \in S : \bar{c}_j \approx v_t\}$ (columns at the threshold)
   4. Let $Q = \{j : \bar{c}_j \geq v_t\}$ (columns at or above threshold)
   5. Find row $i^{\star} \in T$ covered by some column in $J$, choosing the row that minimizes the number of new variables added: $|N_{i^{\star}} \setminus (Q \cup W)|$
   6. Find column $j^{\star} \in J$ covering row $i^{\star}$
   7. Add columns $N_{i^{\star}} \setminus Q$ to the cut set $W$
   8. Update dual progress: $y \leftarrow y + \bar{c}_{j^{\star}}$
   9. For all $j \in N_{i^{\star}} \cap Q$: reduce $\bar{c}_j \leftarrow \bar{c}_j - \bar{c}_{j^{\star}}$
   10. Remove $j^{\star}$ from $S$
3. If $y \geq z^*$: return cover cut $\sum_{j \in W} x_j \geq 1$, else no cut found

**Intuition:** The algorithm greedily builds a cover by exploiting the dual structure. At each step, it identifies a tight constraint whose covering columns have high reduced costs, adds the cheapest columns (in the reduced-cost sense) to the cut, and subtracts their contribution from remaining candidates. The procedure succeeds when the accumulated reduced cost reaches the incumbent bound, guaranteeing the cut is valid.

### 4.3 Per-Row Cover Heuristic

When BCG fails to produce a cut, the `BalasCutSeparator` falls back to a simpler heuristic:

1. For each row $i$ with positive dual value $u_i > \epsilon$:
   - Collect $W_i = \{j : A_{ij} > 0, \bar{x}_j > 0, \bar{c}_j < -\epsilon\}$
   - Check violation: $1 - \sum_{j \in W_i} \bar{x}_j > \epsilon$
   - If violated, create cover cut $\sum_{j \in W_i} x_j \geq 1$
2. Sort candidates by violation magnitude (descending)
3. Return up to 50 cuts

This is a fast fallback that catches simple violated covers from individual constraints.

## 5. Balas Branch Generation (BBG)

### 5.1 Multi-Way Branching

Beyond cut generation, the Balas additive procedure also generates _branch sets_ for multi-way branching. The BBG algorithm (`balas_branch_generate()`) extends BCG to produce multiple disjoint sets $R_1, R_2, \ldots, R_p$, each defining a branch child.

**Algorithm (similar to BCG but accumulates sets):**

The iteration loop is analogous to BCG, but instead of building a single cut, each iteration produces a branch set $R_k = N_{i^*} \cap Q$ — the columns covering the selected tight row that are at or above the current threshold. The algorithm continues until the accumulated reduced cost reaches the incumbent bound or all candidates are exhausted.

**BR1 Multi-Way Branching:**

The solver creates $p$ children from $p$ branch sets:

- Child $k$ fixes all variables in $R_k$ to 0 (none of these columns selected)
- Child $k$ adds cover constraints for all earlier sets: $\sum_{j \in R_i} x_j \geq 1$ for $i < k$

This partitions the search space: child $k$ corresponds to the scenario where sets $R_1, \ldots, R_{k-1}$ are all "covered" (at least one column selected from each) but $R_k$ is not.

**Activation Criteria:**

BR1 branching is only used when the sets provide sufficient diversification:

- At least 3 sets produced ($p \geq 3$)
- Total variables across sets exceed $p \log_2 p$ (diversity check)
- At most 1 singleton set (sets of size 1 are just regular binary branching)
- Total variables at least $2p$ (sufficient branching power)
- $p$ within configured limits (`balas_max_branches`, capped at 12 normally, 8 when forced)

### 5.2 Stagnation-Triggered Balas

During branch-and-bound, if the solver detects _stagnation_ (no gap improvement for `gap_stagnation_window` nodes), it triggers aggressive Balas branching:

1. A frequency accumulator increments by `aggressive_balas_frequency` (default 0.6) at each stagnation event
2. When the accumulator reaches 1.0, the next fractional node forces BBG with `force_balas = true`
3. If BR1 criteria are met, the node creates multiple children via the Balas branch sets
4. Cover cuts from the branch sets are also added to strengthen the LP

This mechanism helps escape plateaus where standard binary branching makes little progress.

## 6. Cut Management

### 6.1 Root Cut Rounds

At the root node, the solver performs up to `cut_rounds_root` (default 5) rounds of cut separation:

1. Solve LP relaxation
2. For each separator (dual aggregated CG, row-pair CG, Balas), generate violated cuts
3. Append cuts to the base model
4. Rebuild LP with cuts and resolve
5. Check if dual bound improved; if no cuts found, stop early

After the standard cut rounds, an additional Balas cover cut phase generates branch-derived covers from BBG.

**Cut Validation for Integer Objectives:**

When all objective coefficients are integer, the LP bound can be tightened to $\lceil z_{LP} \rceil$. If the cut rounds improve the raw LP bound but the tightened (ceiled) bound doesn't change, the cuts are _rolled back_ — they add LP rows without actually improving the effective bound, which wastes LP solve time in the BnB tree.

### 6.2 Mid-BnB Cuts

Mid-BnB cut generation is controlled by `mid_bnb_cut_frequency` (default 0.0 = disabled). When enabled, cuts are generated at stagnation events using the same separators as the root, with the same validation and rollback logic.

### 6.3 Cut Storage

Cuts are stored at two levels:

- **Base model cuts** (`base.base_cuts`): persistent across all BnB nodes, added to the LP at every node solve
- **Node cuts** (`node.cuts`): local to a subtree, used for Balas branch cover constraints

When columns are eliminated during mid-BnB preprocessing, all cuts are remapped to the new variable indices. Cuts whose variables are entirely eliminated are removed.

## 7. Configuration

| Parameter                    | Default | Description                                 |
| ---------------------------- | ------- | ------------------------------------------- |
| `cuts_enabled`               | `true`  | Enable cut separation                       |
| `cut_rounds_root`            | `5`     | Maximum cut rounds at root                  |
| `max_cuts_per_round`         | `100`   | Maximum cuts added per round                |
| `mid_bnb_cut_frequency`      | `0.0`   | Mid-BnB cut frequency (0 = disabled)        |
| `mid_bnb_cut_rounds`         | `3`     | Rounds of mid-BnB cutting                   |
| `balas_enabled`              | `true`  | Enable Balas branching                      |
| `aggressive_balas_frequency` | `0.6`   | Stagnation Balas trigger frequency          |
| `balas_max_branches`         | `20`    | Maximum BR1 children                        |
| `gap_stagnation_window`      | `50`    | Nodes without improvement before stagnation |

## 8. References

1. E. Balas and A. Ho. Set covering algorithms using cutting planes, heuristics, and subgradient optimization. _Mathematical Programming Study_, 12:37–60, 1980.
2. V. Chvátal. Edmonds polytopes and a hierarchy of combinatorial problems. _Discrete Mathematics_, 4(4):305–337, 1973.
3. R. E. Gomory. Outline of an algorithm for integer solutions to linear programs. _Bulletin of the American Mathematical Society_, 64(5):275–278, 1958.
4. G. L. Nemhauser and L. A. Wolsey. _Integer and Combinatorial Optimization_. Wiley, 1988. Chapter II.2 (Chvátal-Gomory cuts) and Chapter II.3 (cover inequalities).
