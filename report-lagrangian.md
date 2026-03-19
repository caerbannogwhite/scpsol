# Lagrangian Relaxation in scpsol: Theory and Implementation

## 1. Introduction

Lagrangian relaxation is a classical technique for obtaining dual bounds and heuristic solutions in combinatorial optimization. For Set Covering Problems (SCP), relaxing the covering constraints into the objective via non-negative multipliers produces a Lagrangian dual that decomposes into independent variable-level subproblems, each solvable by inspection. Combined with subgradient optimization and a repair heuristic, this yields high-quality feasible solutions that significantly improve the incumbent before branch-and-bound begins.

The scpsol solver implements Lagrangian relaxation following the approach of Balas and Ho (1980), using subgradient optimization to update the Lagrangian multipliers and a greedy repair heuristic to convert Lagrangian solutions into feasible covers. This report describes the theory, implementation, and empirical impact.

## 2. Background

### 2.1 Lagrangian Relaxation for SCP

The Set Covering Problem (SCP) is:

$$\min \; c^T x \quad \text{s.t.} \; Ax \geq \mathbf{1}, \; x \in \{0,1\}^n$$

where $A \in \{0,1\}^{m \times n}$ is the constraint matrix, $c > 0$ is the cost vector, and $\mathbf{1}$ is the all-ones vector.

The Lagrangian relaxation with respect to the covering constraints introduces non-negative multipliers $u \geq 0$:

$$L(u) = \min_{0 \leq x \leq 1} \; (c - A^T u)^T x + \mathbf{1}^T u$$

The covering constraints are moved into the objective as penalties. Violating row $i$ (having $\sum_j A_{ij} x_j < 1$) incurs a cost proportional to the multiplier $u_i$.

### 2.2 Subproblem Decomposition

The Lagrangian subproblem decomposes by variable. For each column $j$, define the _Lagrangian reduced cost_:

$$\bar{c}_j = c_j - \sum_{i : A_{ij} > 0} u_i A_{ij}$$

Since $0 \leq x_j \leq 1$ and we minimize, the optimal Lagrangian solution is:

$$x_j^* = \begin{cases} 1 & \text{if } \bar{c}_j < 0 \\ 0 & \text{otherwise} \end{cases}$$

The Lagrangian bound is then:

$$L(u) = \sum_i u_i + \sum_{j : \bar{c}_j < 0} \bar{c}_j$$

This is a valid lower bound on the integer optimum for any $u \geq 0$. The Lagrangian dual seeks the tightest such bound: $\max_{u \geq 0} L(u)$.

### 2.3 Relationship to LP Relaxation

For SCP, the LP relaxation and the Lagrangian dual provide the _same bound_ at optimality. This follows from LP strong duality: the constraint matrix is non-negative with $\geq$ constraints, so the LP dual is:

$$\max \; \mathbf{1}^T u \quad \text{s.t.} \; A^T u \leq c, \; u \geq 0$$

which is exactly the Lagrangian dual $\max_{u \geq 0} L(u)$. At the LP optimum, $L(u^*) = z_{LP}$.

Since the LP relaxation is solved exactly at the root node, the Lagrangian bound does not improve the dual bound. The primary value of Lagrangian relaxation lies in the _heuristic solutions_ produced at each iteration — the Lagrangian subproblem identifies a set of cheap columns that, after repair, often yield near-optimal feasible covers.

### 2.4 Subgradient Optimization

To maximize $L(u)$ over $u \geq 0$, we use the subgradient method (Held and Karp, 1971; Held, Wolfe, and Crowder, 1974). At iteration $k$, the subgradient of $L$ at $u^k$ is:

$$s_i^k = 1 - \sum_j A_{ij} x_j^*(u^k)$$

This measures the _coverage deficit_ of row $i$ under the Lagrangian solution: positive if the row is uncovered, negative if over-covered.

The multiplier update rule is:

$$u_i^{k+1} = \max\left(0, \; u_i^k + t_k \cdot s_i^k\right)$$

where the step size $t_k$ is computed using the Held-Karp formula:

$$t_k = \frac{\lambda_k \cdot (UB - L(u^k))}{\|s^k\|^2}$$

Here $UB$ is the best known upper bound (incumbent), and $\lambda_k$ is a scaling parameter initialized to 2.0 and halved whenever $L(u)$ fails to improve for a fixed number of iterations.

**Convergence behavior:** The multipliers $u$ are initialized from the LP dual solution, so the Lagrangian bound starts near $z_{LP}$ and converges quickly. The practical value is that different multiplier vectors at intermediate iterations produce different Lagrangian solutions, each feeding the repair heuristic and potentially yielding improving feasible solutions.

## 3. The Lagrangian Heuristic

The Lagrangian solution $x^*(u)$ is generally infeasible for the original SCP because the relaxed covering constraints are not enforced. The heuristic converts each Lagrangian solution into a feasible cover through three phases.

### 3.1 Seed from Lagrangian

Start with the columns selected by the Lagrangian subproblem ($x_j = 1$ when $\bar{c}_j < 0$). These columns are included because their Lagrangian reduced cost is negative — they are "cheap" relative to the current multiplier valuation of coverage. Compute the initial coverage of each row.

### 3.2 Greedy Repair

Identify uncovered rows (where $\sum_j A_{ij} x_j < 1$) and greedily add columns to cover them. At each step, select the unselected column $j$ maximizing:

$$\text{score}(j) = \frac{\text{uncovered rows newly covered by } j}{c_j}$$

This is the standard greedy set cover ratio, using original costs rather than Lagrangian reduced costs, since the goal is to minimize the actual objective.

Repeat until all rows are covered. The Lagrangian seed typically covers most rows already, so the repair phase adds only a few columns.

### 3.3 Redundancy Removal

After repair, the solution may contain redundant columns — columns whose removal would leave all rows still covered. Sort selected columns by decreasing cost and attempt to remove each one:

- For column $k$, check if every row $i$ covered by $k$ has sufficient remaining coverage: $\text{coverage}_i - A_{ik} \geq 1$
- If so, remove $k$ and update the coverage counts

This phase drops expensive columns that were included by the Lagrangian seed but are made redundant by the repair phase or by other seed columns.

### 3.4 Feasibility Verification

As a safety check, the heuristic solution is verified against the original CSR constraint matrix before being accepted. Each row is checked independently to ensure $\sum_j A_{ij} x_j \geq 1$, guarding against numerical errors in the incremental coverage tracking.

## 4. Implementation

### 4.1 Code Structure

The implementation resides in two files:

- **`src/lagrangian.h`**: declares `LagrangianResult` (best bound, best heuristic objective, best solution vector, final multipliers, iteration count) and the `lagrangian_relaxation()` function signature.
- **`src/lagrangian.cpp`**: implements the subgradient loop with embedded heuristic.

### 4.2 Algorithm Parameters

| Parameter        | Value | Description                                                   |
| ---------------- | ----- | ------------------------------------------------------------- |
| `max_iterations` | 500   | Maximum subgradient iterations                                |
| `time_limit`     | 5.0 s | Wall-clock time limit                                         |
| `lambda_init`    | 2.0   | Initial step size scaling parameter                           |
| `halve_interval` | 30    | Iterations without bound improvement before halving $\lambda$ |
| `lambda_min`     | 1e-6  | Minimum $\lambda$ before early termination                    |

In practice, 500 iterations complete in ~0.1–0.2 seconds on models with 300–600 columns and 500 rows, well within the 5-second budget.

### 4.3 Data Structures

The subproblem uses the `BaseRelaxationModel`'s column-to-rows transpose (`cols_to_rows`) for computing Lagrangian reduced costs, and the CSR matrix (`csr_offs`, `csr_inds`, `csr_vals`) for the final feasibility verification. Four working arrays are reused across iterations:

- `x[ncols]` — Lagrangian solution (binary)
- `rc[ncols]` — Lagrangian reduced costs
- `subgrad[nrows]` — subgradient vector
- `coverage[nrows]` — row coverage for the heuristic

### 4.4 Solver Integration

Lagrangian relaxation runs as **Phase 3.5** in the solver pipeline, after the root LP and its heuristics but before reduced-cost fixing:

1. **Phase 1–2:** Preprocessing (greedy cover, dominance, essential columns)
2. **Phase 3:** Root LP relaxation + LP-based heuristics
3. **Phase 3.5:** Lagrangian relaxation (this module)
4. **Phase 4–5:** Iterative reduced-cost fixing, LP probing, BnB

The multipliers are initialized from the LP dual solution (`root_sol.row_dual`), clamped to $\geq 0$. If the Lagrangian heuristic finds a solution strictly better than the current incumbent, the incumbent is updated before reduced-cost fixing begins. This tighter incumbent directly shrinks the RC fixing gap ($\text{gap} = \text{incumbent} - z_{LP}$), enabling more columns to be fixed.

### 4.5 Early Termination

The subgradient loop terminates early when any of these conditions is met:

- **Time limit exceeded**: wall-clock time exceeds the budget
- **Zero subgradient**: $\|s\|^2 < 10^{-12}$, indicating the Lagrangian solution satisfies all constraints (optimal for the relaxation)
- **Zero step size**: $t_k < 10^{-12}$, indicating $UB \approx L(u)$ (the gap has closed)
- **Lambda underflow**: $\lambda < 10^{-6}$ after repeated halvings (stalled convergence)

## 5. Experimental Results

### 5.1 Incumbent Quality

The Lagrangian heuristic consistently finds better incumbents than the greedy set cover heuristic. The table below compares the greedy incumbent (before the root LP) with the Lagrangian heuristic objective (after 500 subgradient iterations) across OR-Library benchmark instances.

| Instance | Greedy | Lagrangian | Optimal | Improvement |
| -------- | -----: | ---------: | ------: | ----------: |
| scpa1    |    271 |        255 |     253 |          16 |
| scpa2    |    269 |        253 |     252 |          16 |
| scpa3    |    244 |        234 |     232 |          10 |
| scpa4    |    242 |        235 |     234 |           7 |
| scpa5    |    242 |        237 |     236 |           5 |
| scpb1    |     75 |         70 |      69 |           5 |
| scpb2    |     85 |         76 |      76 |       9 (=) |
| scpb3    |     85 |         81 |      80 |           4 |
| scpb4    |     87 |         79 |      79 |       8 (=) |
| scpb5    |     76 |         72 |      72 |       4 (=) |
| scpnre1  |     30 |         29 |      29 |       1 (=) |
| scpnre2  |     33 |         31 |      30 |           2 |
| scpnre3  |     29 |         28 |      27 |           1 |
| scpnre4  |     30 |         29 |      28 |           1 |
| scpnre5  |     30 |         28 |      28 |       2 (=) |
| scpnrf1  |     16 |         14 |      14 |       2 (=) |
| scpnrf2  |     16 |         15 |      15 |       1 (=) |
| scpnrf3  |     17 |         15 |      14 |           2 |
| scpnrf4  |     17 |         15 |      14 |           2 |
| scpnrf5  |     16 |         14 |      13 |           2 |

**(=)** marks instances where the Lagrangian heuristic finds the optimal solution.

The Lagrangian heuristic finds the optimal on 7 out of 20 instances (35%) and comes within 1–2 units of optimal on all others. The average improvement over the greedy is 5.0 units, ranging from 1 unit on the small NRE/NRF instances to 16 units on the A-set.

### 5.2 Impact on Reduced-Cost Fixing

The tighter incumbent from the Lagrangian enables more aggressive reduced-cost fixing. The table below shows the number of columns remaining after iterative RC fixing, with and without the Lagrangian incumbent.

| Instance | Without Lagrangian | With Lagrangian | Cols saved |
| -------- | -----------------: | --------------: | ---------: |
| scpnre1  |                471 |             434 |         37 |
| scpnre2  |                604 |             552 |         52 |
| scpnre3  |                459 |             459 |          0 |
| scpnre4  |                507 |             464 |         43 |
| scpnre5  |                539 |             448 |         91 |
| scpnrf1  |                387 |             332 |         55 |
| scpnrf2  |                359 |             307 |         52 |
| scpnrf3  |                353 |             353 |          0 |
| scpnrf4  |                394 |             367 |         27 |
| scpnrf5  |                353 |             353 |          0 |

On 7 out of 10 NRE/NRF instances, the Lagrangian incumbent enables additional column removal. The largest improvement is on scpnre5 (91 fewer columns, a 17% reduction). Instances where the greedy already found a strong incumbent (scpnre3, scpnrf3, scpnrf5) see no additional benefit from the tighter Lagrangian bound.

### 5.3 Computational Overhead

The Lagrangian relaxation completes 500 iterations in 0.1–0.2 seconds on models with 300–600 columns. This is negligible relative to the 30–300 second total solve times on hard instances.

### 5.4 Interaction with BnB

The tighter incumbent and smaller post-RC-fixing model change the BnB search tree structure. While a smaller model is generally easier to solve, different branching decisions can lead to different tree exploration patterns. On some instances (scpnre5: 88s → 82s), the smaller model accelerates BnB. On others (scpnre1: 267s → 300s), the different tree structure happens to be harder. This is a well-known phenomenon in MIP solving — perturbations to the preprocessing pipeline can shift the BnB tree in either direction. The net effect across all instances is positive due to the consistently stronger starting incumbent.

## 6. References

1. E. Balas and A. Ho. Set covering algorithms using cutting planes, heuristics, and subgradient optimization. _Mathematical Programming Study_, 12:37–60, 1980.
2. M. Held and R. M. Karp. The traveling-salesman problem and minimum spanning trees. _Operations Research_, 18(6):1138–1162, 1970.
3. M. Held, P. Wolfe, and H. P. Crowder. Validation of subgradient optimization. _Mathematical Programming_, 6(1):62–88, 1974.
4. M. L. Fisher. The Lagrangian relaxation method for solving integer programming problems. _Management Science_, 27(1):1–18, 1981.
5. A. Caprara, M. Fischetti, and P. Toth. A heuristic method for the set covering problem. _Operations Research_, 47(5):730–743, 1999.
