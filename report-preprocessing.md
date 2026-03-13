# Preprocessing Techniques in scpsol: Theory and Implementation

## 1. Introduction

Preprocessing (or presolving) reduces the size of a combinatorial optimization problem before the main solver begins. For Set Covering Problems (SCP), this means removing redundant rows (constraints), dominated columns (variables), and fixing variables whose optimal values can be deduced. Effective preprocessing can dramatically reduce the search space and solve time, sometimes resolving the problem entirely without branch-and-bound.

The scpsol solver implements a multi-stage preprocessing pipeline that iterates multiple techniques until a fixpoint is reached — that is, until no further reductions are possible. This report describes the theory behind each technique, the implementation in scpsol, and how they interact.

## 2. Problem Formulation

The Set Covering Problem is:

$$\min \sum_{j=1}^{n} c_j x_j$$
$$\text{subject to } \sum_{j \in S_i} x_j \geq 1, \quad i = 1, \ldots, m$$
$$x_j \in \{0, 1\}, \quad j = 1, \ldots, n$$

where $c_j > 0$ is the cost of column $j$, $S_i$ is the set of columns that cover row $i$, and the constraint requires every row to be covered by at least one selected column. We denote the set of rows covered by column $j$ as $R_j$.

## 3. Preprocessing Techniques

### 3.1 Essential Column Fixing

**Theory.** A column $j$ is *essential* for row $i$ if it is the only column that can cover $i$, i.e., $S_i = \{j\}$. In any feasible solution, $x_j = 1$ must hold. Fixing $x_j = 1$ allows removing all rows in $R_j$ from the problem, since they are now guaranteed to be covered.

Crucially, removing rows can make other columns essential. If column $k$ was one of two columns covering row $i'$, and the other column was just fixed, then $k$ becomes essential for $i'$. This creates a *cascade* of fixings.

**Implementation.** The essential column cascade in `row_reduce()` (`src/preprocessor.cpp`) iterates until quiescence:

```
repeat:
    for each active row i:
        count active columns covering i
        if count == 1:
            fix that column to 1, remove all its rows
            mark changed = true
until not changed
```

The fixed columns' costs are accumulated in `fixed_preprocess_cost` and subtracted from the incumbent bound, since they are known to be in every optimal solution.

### 3.2 Row Domination

**Theory.** Row $i$ *dominates* row $i'$ if every column covering $i$ also covers $i'$, i.e., $S_i \subseteq S_{i'}$. Any feasible solution that satisfies the constraint for $i$ automatically satisfies the constraint for $i'$, so $i'$ is redundant and can be removed.

**Implementation.** In `row_reduce()`, active rows are sorted by coverage size (ascending), and for each pair $(i, i')$ where $|S_i| \leq |S_{i'}|$, a sorted-set subset check determines if $S_i \subseteq S_{i'}$. If so, $i'$ is deactivated.

### 3.3 Probing (Constraint Propagation)

**Theory.** For each active column $j$, we temporarily assume $x_j = 0$ and propagate the consequences via the essential column cascade. If this leads to a row with no remaining covering columns (infeasibility), then $x_j = 0$ is impossible in any feasible solution, so $x_j = 1$ must hold. This is equivalent to a single round of *node presolve* or *probing* in MIP solvers.

More formally: setting $x_j = 0$ may force other columns to be essential (creating $x_k = 1$), which removes rows and may make further columns essential. If the cascade reaches a contradiction (an uncoverable row), the original assumption $x_j = 0$ is refuted.

**Implementation.** In `row_reduce()`, for each active column $j$:

1. Temporarily deactivate $j$
2. Run the essential column cascade, tracking all changes on a rollback stack
3. If any row becomes uncoverable, mark $j$ as fixed to 1
4. Rollback all temporary changes

After all probing, a final essential column cascade picks up any newly created essentials.

### 3.4 Single Column Dominance

**Theory.** Column $j$ is *dominated* by column $k$ if $c_k \leq c_j$ and $R_j \subseteq R_k$ — that is, $k$ covers all the rows that $j$ covers, at no greater cost. In any optimal solution that uses $j$, replacing $j$ with $k$ yields a solution that is at least as good. Therefore $j$ can be removed.

**Implementation.** The `SingleColumnDominanceRule` checks every pair of active columns. For each target column $j$, if any candidate $k$ has $c_k \leq c_j$ and $R_j \subseteq R_k$ (verified via sorted-set inclusion), then $j$ is deactivated. Tie-breaking: if costs are equal, the lower-indexed column is kept.

### 3.5 Two-Column Dominance

**Theory.** Column $j$ is *pair-dominated* if there exist columns $k_1, k_2$ such that $c_{k_1} + c_{k_2} < c_j$ and $R_j \subseteq R_{k_1} \cup R_{k_2}$. Any solution using $j$ can be improved by using $k_1$ and $k_2$ instead, so $j$ is redundant.

**Implementation.** The `TwoColumnDominanceRule` iterates over all target columns and all pairs of candidates. The `unionCoversSorted()` helper efficiently checks whether the union of two sorted row sets covers the target's rows. Time limit checks prevent excessive computation on large instances.

### 3.6 Cost-Driven Replacement

**Theory.** This generalizes pair dominance by considering 2-column and 3-column replacements, focused on expensive columns first. By processing columns in descending cost order, the most impactful reductions are found first. A column $j$ is removed if any combination of 2 or 3 cheaper columns covers all its rows at lower total cost.

**Implementation.** The `CostDrivenReplacementRule` builds a reverse index (columns by row) to efficiently find candidate replacements. For each target column (processed from most expensive to cheapest), it collects all columns sharing at least one row with the target, sorts them by cost, and checks pairs and triples.

### 3.7 Incumbent Budget Pruning

**Theory.** Given a known upper bound $z^*$ (incumbent), any column $j$ with $c_j \geq z^*$ cannot appear in any improving solution, since it alone would exceed the budget. More generally, if including $j$ requires a minimum additional cost to cover the remaining rows that exceeds $z^* - 1$, then $j$ can be pruned.

The budget for column $j$ is $B_j = \lfloor z^* \rfloor - 1 - \lfloor c_j \rfloor$. If $B_j < 0$, the column is immediately pruned. For small budgets ($B_j < 1$), we check whether a single cost-1 column can cover all rows not covered by $j$. For larger budgets, we compute a lower bound on the cost of covering the remaining rows and prune if it exceeds $B_j$.

**Implementation.** The `IncumbentBudgetPruningRule` processes columns from most expensive to cheapest, building a per-row sorted column list for efficient minimum-cost lookups.

### 3.8 Reduced-Cost Fixing

**Theory.** After solving the LP relaxation with dual bound $z_{LP}$ and incumbent $z^*$, the *gap* is $\Delta = z^* - z_{LP}$. For any column $j$ with LP value $x_j = 0$ and reduced cost $\bar{c}_j > \Delta$, fixing $x_j = 1$ would increase the LP bound beyond $z^*$, so $x_j = 0$ in every optimal solution. This is a standard technique in mixed-integer programming.

**Implementation.** In `solver.cpp` (Phase 4-5), reduced costs are computed from the LP dual solution. Columns at their lower bound ($x_j \approx 0$) with $\bar{c}_j > \Delta + \epsilon$ are removed.

### 3.9 Cost Reduction

**Theory.** Any column with $c_j \geq z^*$ (the incumbent bound) can be removed, since including it alone would match or exceed the best known solution. Combined with the SCP constraint that at least one other column must be present, the total would strictly exceed $z^*$.

**Implementation.** The `apply_cost_reduction()` function removes all columns with cost at or above the incumbent bound.

## 4. Greedy Heuristic

### 4.1 Theory

The greedy algorithm by Chvátal (1979) provides an initial feasible solution and upper bound. At each iteration, it selects the column with the best *cost-effectiveness ratio*:

$$j^* = \arg\max_{j} \frac{|\{i \in R_j : i \text{ uncovered}\}|}{c_j}$$

That is, it picks the column that covers the most uncovered rows per unit cost. This has a well-known approximation guarantee: the greedy solution costs at most $H(m)$ times the optimal, where $H(m) = 1 + 1/2 + \cdots + 1/m$ is the $m$-th harmonic number.

### 4.2 Implementation

The `greedy_set_cover_heuristic()` function in `src/preprocessor.cpp` implements the Chvátal greedy:

1. Build a reverse index (rows covered by each column)
2. Maintain a set of uncovered rows
3. At each step, scan all unused columns, compute the number of uncovered rows each covers, and select the one with the highest coverage/cost ratio
4. Mark covered rows, add the column to the solution
5. Repeat until all rows are covered

The greedy solution provides the initial incumbent bound $z^*$ that enables cost reduction and budget pruning in subsequent preprocessing stages.

## 5. Preprocessing Pipeline

### 5.1 Iterated Fixpoint

A key design principle in scpsol is that preprocessing techniques are *iterated until fixpoint*. Each technique can enable further reductions by others:

- **Essential columns → dominance**: Fixing essential columns and removing their rows can make previously non-dominated columns now dominated.
- **Dominance → essential columns**: Removing dominated columns can leave rows with only one covering column, creating new essentials.
- **Probing → dominance**: Probing may fix columns and remove rows, enabling new dominance relationships.
- **Budget pruning → essential columns**: Removing expensive columns can create new essential columns.

The preprocessing loop runs all techniques in sequence and repeats the entire sequence until no further reductions occur (or a maximum of 10 rounds).

### 5.2 Pipeline Stages

The full preprocessing pipeline in `solver.cpp` proceeds as follows:

**Phase 1: Greedy Heuristic**
- Compute initial feasible solution and upper bound $z^*$

**Phase 2: Pre-LP Reduction (iterated)**
Repeat until fixpoint:
1. Cost reduction (remove columns with $c_j \geq z^*$)
2. Incumbent budget pruning
3. Cost-driven replacement (2- and 3-column dominance)
4. Configured dominance rules (single, two-column)
5. Row reduction (essential columns, row domination, probing)

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

## 6. Configuration

The preprocessing behavior is controlled by the `SolverConfig` fields:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `preprocess_rules` | `"single,two"` | Comma-separated list of dominance rules |
| `preprocess_time_limit` | `10.0` | Time limit (seconds) per preprocessing phase |

Available rule tokens: `single` (single column dominance), `two` (two-column dominance), `cost_driven` (cost-driven replacement), `incumbent_budget` (budget pruning), `none` (disable all).

## 7. References

- V. Chvátal. A greedy heuristic for the set-covering problem. *Mathematics of Operations Research*, 4(3):233–235, 1979.
- E. Balas and A. Ho. Set covering algorithms using cutting planes, heuristics, and subgradient optimization. *Mathematical Programming Study*, 12:37–60, 1980.
- T. Achterberg, T. Koch, and A. Martin. Constraint integer programming: a new approach to integrate CP and MIP. *LNCS*, 3011:6–20, 2004.
- M. Savelsbergh. Preprocessing and probing techniques for mixed integer programming problems. *ORSA Journal on Computing*, 6(4):445–454, 1994.
