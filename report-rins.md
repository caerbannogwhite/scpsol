# RINS Heuristic in scpsol: Theory and Implementation

## 1. Introduction

Relaxation Induced Neighborhood Search (RINS) is a primal improvement heuristic for mixed-integer programs introduced by Danna, Rothberg, and Le Pape (2005). The key insight is that variables where the LP relaxation and the incumbent integer solution agree are likely to retain their values in the optimal solution. By fixing these "consensus" variables and solving the residual problem on the remaining "disagreement" variables, RINS explores a focused neighborhood around the incumbent that often contains better solutions.

The scpsol solver implements RINS as a periodic heuristic during branch-and-bound, parameterized by a frequency accumulator. On the OR-Library benchmark, RINS at low frequencies (0.005–0.01) finds improved incumbents on hard NRF instances that no other heuristic discovers, including the optimal solution for scpnrf5.

## 2. Background

### 2.1 Motivation

The Set Covering Problem (SCP) is $\min \; c^T x$ subject to $Ax \geq \mathbf{1}$, $x \in \{0,1\}^n$, where $A \in \{0,1\}^{m \times n}$ is the constraint matrix and $c > 0$ is the cost vector.

During branch-and-bound, the solver maintains two key solutions:

- **LP relaxation** $\bar{x}$: fractional, optimal for the relaxed problem at the current node
- **Incumbent** $x^*$: the best known integer-feasible solution

These two solutions often agree on most variables: if $\bar{x}_j \approx 1$ and $x^*_j = 1$, column $j$ is very likely part of the optimal cover. If $\bar{x}_j \approx 0$ and $x^*_j = 0$, column $j$ is very likely excluded. The variables where they disagree — where the LP "wants" a different value than the incumbent — define a small neighborhood that concentrates the potential for improvement.

### 2.2 Neighborhood Definition

Given the LP solution $\bar{x}$ and incumbent $x^*$, RINS classifies each variable:

- **Fixed to 1:** $\bar{x}_j > 1 - \delta$ and $x^*_j = 1$ (both agree on inclusion)
- **Fixed to 0:** $\bar{x}_j < \delta$ and $x^*_j = 0$ (both agree on exclusion)
- **Free:** all other variables (LP and incumbent disagree)

where $\delta$ is a tolerance parameter (0.1 in scpsol). The free variables define the RINS neighborhood.

For SCP, this decomposition has a natural interpretation: the "fixed to 1" columns form a partial cover that both the LP and incumbent agree on. The remaining uncovered rows must be covered by the free columns, forming a smaller residual SCP.

## 3. Algorithm

### 3.1 RINS Procedure

Given the current LP solution, incumbent, and base model:

1. **Classify variables:** For each column $j$, determine if it is fixed to 1, fixed to 0, or free based on LP-incumbent agreement.

2. **Guard conditions:** Skip RINS if:
   - Fewer than $n/3$ variables are fixed (insufficient agreement — the sub-problem would be nearly as large as the original)
   - More than $n/2$ variables are free (sub-problem too large)
   - No incumbent exists

3. **Identify covered rows:** Rows covered by any fixed-to-1 column are already satisfied and excluded from the sub-problem.

4. **Build residual SCP:** Construct a sub-problem containing only:
   - Rows not covered by fixed-to-1 columns
   - Free columns that cover at least one uncovered row

5. **Solve sub-problem:** Run the full solver on the residual SCP with reduced limits (5 seconds, 5000 nodes).

6. **Combine solutions:** The RINS solution is the union of fixed-to-1 columns and the sub-problem's optimal cover. If the combined cost improves the incumbent, update it.

### 3.2 Sub-Problem Construction

The residual SCP is packaged as an independent `ScpInstance` with:

- **Rows:** uncovered rows (not covered by any consensus-included column)
- **Columns:** free columns that are useful (cover at least one uncovered row)
- **Costs:** inherited from the parent model
- **Index mapping:** maintained for solution reconstruction

The sub-problem typically has 10–50% of the original variables and 20–80% of the rows, depending on the degree of LP-incumbent agreement. After preprocessing within the sub-solver, it is often much smaller.

### 3.3 Budget Control

RINS calls the full `solve()` function on the sub-problem, which includes preprocessing, root LP, Lagrangian heuristic, and branch-and-bound. To bound the overhead:

- **Time limit:** 5 seconds per RINS invocation
- **Node limit:** 5000 nodes
- **No recursion:** RINS, diving, and decomposition are disabled in the sub-solver
- **Frequency accumulator:** RINS triggers every $\lfloor 1/f \rfloor$ nodes, where $f$ is the configured frequency

At frequency 0.01, RINS fires approximately once every 100 nodes. With a 5-second budget per call, the total RINS overhead on a 20000-node solve is at most $200 \times 5 = 1000$ seconds — but in practice, most sub-problems solve in well under 1 second due to their small size.

## 4. Implementation

### 4.1 Integration with BnB

RINS is called after the existing heuristic block and before pruning-by-bound:

1. LP solve produces $\bar{x}$ and duals
2. Standard heuristics (nearest integer, dual-guided repair)
3. Lagrangian heuristic (if enabled)
4. Diving heuristic (if enabled)
5. **RINS** (if enabled and not in proving phase)
6. Pruning by bound

The incumbent in active column space is reconstructed lazily before each RINS call by reverse-mapping through `active_to_original`.

### 4.2 Proving Phase Suppression

RINS is a primal improvement heuristic — it searches for better incumbents. When the BnB tree is in the "proving phase" (frontier monotonically shrinking, meaning the tree is closing rather than expanding), RINS cannot help because the incumbent is already optimal or near-optimal. The solver detects this phase via a hysteresis-based frontier trend tracker:

- **Enter proving phase:** after 5 consecutive stagnation windows where the frontier shrank
- **Exit proving phase:** after 3 consecutive stagnation windows where the frontier grew

This avoids the oscillation that occurs with a tighter window. During proving phase, RINS, diving, and mid-BnB cuts are all suppressed, allowing the solver to focus LP budget on node processing.

### 4.3 Frequency Selection

The frequency parameter operates on a logarithmic scale in practice:

| Frequency | Approx. period | Overhead    | Use case                              |
| --------- | -------------- | ----------- | ------------------------------------- |
| 0.001     | ~1000 nodes    | Negligible  | Background improvement, minimal cost  |
| 0.005     | ~200 nodes     | Low         | Good balance for hard instances       |
| 0.01      | ~100 nodes     | Moderate    | Aggressive improvement, more overhead |
| 0.1       | ~10 nodes      | High        | Very aggressive, may slow BnB         |

The optimal frequency depends on instance hardness: for instances that solve quickly (< 100 nodes), any RINS frequency adds pure overhead. For hard instances (> 10000 nodes) where the incumbent is suboptimal, frequencies of 0.005–0.01 provide the best tradeoff.

## 5. Experimental Results

### 5.1 Setup

All experiments use the OR-Library SCP benchmark instances (Beasley, 1990) with a 300-second time limit, single-threaded, on Windows 11. The baseline is v0.10 with all default settings (RINS disabled).

### 5.2 Results on Core Instances (scpa, scpb)

RINS has no impact on the easy instances: scpa1–5 and scpb1–5 all solve optimally with identical node counts and times across all frequencies. The existing heuristics (greedy + Lagrangian + dual-guided repair) already find the optimal incumbent within the first few nodes.

### 5.3 Results on Hard Instances (scpnre, scpnrf)

| Instance    | Baseline       | RINS 0.001     | RINS 0.005     | RINS 0.01          |
| ----------- | -------------- | -------------- | -------------- | ------------------ |
| scpnre1     | 29 Opt 245s    | 29 Opt 247s    | 29 Opt 265s    | 29 Opt 278s        |
| scpnre2     | 30/29 TL       | 30/29 TL       | 30/29 TL       | 30/29 TL           |
| scpnre3     | 27 Opt 111s    | 27 Opt 107s    | 27 Opt 114s    | 27 Opt 122s        |
| scpnre4     | 28 Opt 219s    | 28 Opt 220s    | 28 Opt 241s    | 28 Opt 254s        |
| scpnre5     | 28 Opt 80s     | 28 Opt 80s     | 28 Opt 83s     | 28 Opt 85s         |
| scpnrf1     | 14 Opt 155s    | 14 Opt 155s    | 14 Opt 162s    | 14 Opt 167s        |
| scpnrf2     | 15 Opt 117s    | 15 Opt 117s    | 15 Opt 123s    | 15 Opt 127s        |
| scpnrf3     | 14 Opt 141s    | 14 Opt 142s    | —              | **14 Opt 34s**     |
| scpnrf4     | 14/13 TL       | 14/13 TL       | 14/13 TL       | 14/13 TL           |
| **scpnrf5** | **14**/13 TL   | **14**/12 TL   | **13**/12 TL   | **13**/12 TL       |

### 5.4 Key Findings

**scpnrf5 improvement.** The most significant result: RINS at frequencies 0.005 and 0.01 finds the optimal incumbent of **13** on scpnrf5, where the baseline (and all other heuristics including greedy, Lagrangian, dual-guided, local search, and diving) find only 14. The LP relaxation at certain BnB nodes provides enough structural guidance for RINS to discover the one-column swap that improves the cover.

**scpnrf3 speedup.** At frequency 0.01, scpnrf3 solves in **34 seconds** with only **2297 nodes**, compared to 141 seconds and 14870 nodes without RINS. The RINS-found incumbent enables aggressive RC fixing and budget pruning that dramatically reduces the remaining tree.

**Overhead scaling.** At 0.001, LP solve overhead is 5–10% (negligible). At 0.01, overhead is 30–50% on NRE instances (moderate), but this is offset by tree reduction when RINS finds improvements. At 0.1, overhead exceeds 100% and the solver processes far fewer nodes within the time limit.

**NRE instances unaffected.** On scpnre1–5, existing heuristics already find optimal incumbents early. RINS searches the LP-incumbent neighborhood but finds no improvements. The overhead from 0.005–0.01 adds 5–15% to solve time without benefit.

## 6. Relationship to Other Techniques

### 6.1 Local Search

scpsol's local search heuristic performs column swaps (remove one column, add a cheaper one that covers the same rows). RINS explores a fundamentally different neighborhood: it re-optimizes a subset of variables rather than performing local moves. The two are complementary — local search finds nearby solutions via single swaps, while RINS can discover multi-column rearrangements.

### 6.2 Diving Heuristics

Diving fixes variables one at a time, re-solving LPs along a single path. RINS fixes many variables at once (all consensus variables) and solves the residual problem optimally. Diving is cheaper per call but explores a narrower neighborhood. RINS is more expensive but can find solutions that diving misses.

### 6.3 Lagrangian Heuristic

The Lagrangian heuristic operates in the dual space, constructing solutions from subgradient-guided column selection. RINS operates in the primal space, using LP-incumbent agreement to define a search neighborhood. On scpnrf5, the Lagrangian heuristic finds primal=14 (matching the greedy incumbent) while RINS finds primal=13 — demonstrating that the two approaches access different solution neighborhoods.

## 7. References

1. E. Danna, E. Rothberg, and C. Le Pape. Exploring relaxation induced neighborhoods to improve MIP solutions. _Mathematical Programming_, 102(1):71–90, 2005.
2. E. Balas and A. Ho. Set covering algorithms using cutting planes, heuristics, and subgradient optimization. _Mathematical Programming Study_, 12:37–60, 1980.
3. T. Achterberg and T. Berthold. Improving the feasibility pump. _Discrete Optimization_, 4(1):77–86, 2007.
4. M. Fischetti and A. Lodi. Local branching. _Mathematical Programming_, 98(1–3):23–47, 2003.
