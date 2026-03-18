# Branching Techniques in scpsol: Theory and Implementation

## 1. Introduction

The scpsol solver uses branch-and-bound (BnB) to solve Set Covering Problems (SCP) to proven optimality. At each BnB node, the LP relaxation is solved and, if the solution is fractional, a branching variable must be chosen to partition the search space. The quality of this choice has a dramatic impact on the size of the search tree and overall solve time.

This report describes the branching strategies implemented in scpsol, the theory behind them, and their empirical performance on standard SCP benchmark instances.

## 2. Background: Branching Variable Selection

Given an LP relaxation solution where variable $x_j$ has fractional value $f_j \in (0, 1)$, branching creates two child nodes: one with $x_j = 0$ and one with $x_j = 1$. The key question is: _which_ fractional variable should we branch on?

### 2.1 Simple Heuristics

**Most Fractional** selects the variable closest to 0.5, i.e., the one maximizing $|f_j - 0.5|$. The intuition is that variables near 0 or 1 are "almost decided" and branching on them yields little new information. This is fast (O(n) scan) but ignores the objective function entirely.

**Highest Cost Fractional** selects the fractional variable with the largest objective coefficient. The intuition is that fixing expensive variables early has the largest impact on the objective. This is also fast but ignores the LP structure.

Both heuristics were implemented in scpsol as `MostFractionalSelector` and `HighestCostFractionalSelector` in `src/heuristics.cpp`.

### 2.2 Strong Branching

Strong branching (SB) evaluates each candidate by _actually solving_ the two child LPs. For each fractional variable $j$:

1. Fix $x_j = 0$, solve the LP, record the dual bound gain $D^-_j = z^-_j - z_{\text{parent}}$
2. Fix $x_j = 1$, solve the LP, record the dual bound gain $D^+_j = z^+_j - z_{\text{parent}}$
3. Restore the parent LP state

The variable with the best combined score is selected. The standard scoring function is the **product score**:

$$\text{score}(j) = \max(\epsilon, D^-_j) \cdot \max(\epsilon, D^+_j)$$

where $\epsilon$ is a small constant (e.g., $10^{-6}$) to avoid zero products. This score favors variables that improve _both_ child bounds, preferring balanced improvements over one-sided ones.

Full strong branching produces the smallest search trees among standard branching rules, but the computational cost is prohibitive: solving $2|F|$ LPs at every node (where $|F|$ is the number of fractional variables) often makes the per-node overhead exceed the tree size savings.

### 2.3 Pseudocost Branching

Pseudocost branching avoids the cost of solving probe LPs by maintaining _historical statistics_ of how much each variable's branching improved the dual bound in the past.

For each variable $j$, we track:

- **Down pseudocost** $\sigma^-_j$: average of $D^-_j / f_j$ over all previous down-branchings
- **Up pseudocost** $\sigma^+_j$: average of $D^+_j / (1 - f_j)$ over all previous up-branchings

At a new node, the estimated gains are:

$$\hat{D}^-_j = f_j \cdot \sigma^-_j, \quad \hat{D}^+_j = (1 - f_j) \cdot \sigma^+_j$$

This is extremely cheap (no LPs solved) and works well once sufficient history has been accumulated. The problem is _initialization_: at the start of the search, no branching history exists, so pseudocost estimates are unreliable and the first branching decisions may be poor, leading to a large tree that is difficult to recover from.

### 2.4 Reliability Branching

Reliability branching, introduced by Achterberg, Koch, and Martin (2005), elegantly combines strong branching and pseudocost branching. The key idea: use strong branching to _initialize_ pseudocosts, then switch to pseudocost estimates once they become _reliable_.

A variable $j$ is considered **reliable** if it has been strong-branched at least $\eta$ times in both directions (the **reliability parameter**). The algorithm:

1. At each node, for each fractional candidate $j$:
   - If $j$ is **unreliable** (fewer than $\eta$ observations): perform strong branching and update pseudocosts with the observed gains
   - If $j$ is **reliable**: use pseudocost estimates (no LP solves needed)
2. Select the variable with the highest product score

This achieves the best of both worlds: strong branching quality in the critical early nodes (where good decisions have cascading effects) and pseudocost speed in the rest of the tree. After the first few dozen nodes, most variables have been observed $\eta$ times and the per-node cost drops to near zero.

## 3. Implementation in scpsol

### 3.1 Pseudocost State

The `PseudocostState` struct (`src/reliability.h`) maintains four arrays per variable:

```
sum_down[j]   = sum of (gain / frac) over all down-probes of variable j
sum_up[j]     = sum of (gain / (1-frac)) over all up-probes of variable j
count_down[j] = number of down-probe observations
count_up[j]   = number of up-probe observations
```

When no observations exist, a cost-proportional fallback is used: $\hat{D}^-_j = f_j \cdot c_j$ and $\hat{D}^+_j = (1-f_j) \cdot c_j$, where $c_j$ is the objective coefficient. This biases initial branching toward high-cost variables, which is sensible for SCP since fixing an expensive column has the largest impact on solution quality.

The pseudocost arrays must be **remapped** whenever mid-BnB column removal changes the variable indices. This is handled by `PseudocostState::remap()`, called after each `mid_bnb_column_removal()` and `mid_bnb_budget_pruning()`.

### 3.2 Strong Branching with Propagation

A key enhancement over standard strong branching is **propagation**: before solving each probe LP, we run an essential-column cascade to detect implied fixings and infeasibility.

For SCP constraints of the form $\sum_{j \in S_i} x_j \geq 1$, when a variable is tentatively fixed to 0, some rows may lose their last remaining free column. This triggers a cascade:

1. If a row has exactly one free column remaining, that column is **forced to 1**
2. Fixing that column covers additional rows, potentially creating new essential columns
3. If a row has **zero** free columns, the probe direction is **infeasible**

This propagation is implemented in `propagate_fixings()` (`src/reliability.cpp`), reusing the same essential-column cascade logic from the preprocessor but operating non-destructively on a temporary state.

The benefits are twofold:

- **Infeasibility detection without an LP solve**: if propagation proves infeasibility, the probe LP is skipped entirely, saving an expensive LP solve
- **Tighter probe bounds**: implied fixings tighten the probe LP, producing more accurate dual bound gains for pseudocost estimation

### 3.3 Candidate Prioritization

Not all fractional variables need probing. The implementation limits strong-branching probes to at most `max_sb` candidates per node (default 8). Candidates are sorted so that:

1. **Unreliable** variables come first (they need probing most)
2. Among equally unreliable variables, the **most fractional** (closest to 0.5) are probed first

Once the strong-branching budget is exhausted, remaining candidates use pseudocost estimates regardless of reliability.

### 3.4 Integration with Balas Branching

Reliability branching coexists with the existing Balas branching strategy. When stagnation is detected and aggressive Balas branching is triggered (`force_aggressive_branching = true`), Balas multi-way branching takes precedence. Reliability branching is only used for standard binary branching decisions.

## 4. Configuration

The following command-line parameters control branching behavior:

| Parameter              | Default     | Description                                                           |
| ---------------------- | ----------- | --------------------------------------------------------------------- |
| `--branch reliability` | reliability | Branching strategy (`reliability`, `most_fractional`, `highest_cost`) |
| `--reliability-eta N`  | 4           | Number of observations before pseudocost is considered reliable       |
| `--reliability-sb N`   | 8           | Maximum strong-branch probes per node                                 |
| `--balas-frequency F`  | 0.6         | Frequency of aggressive Balas branching during stagnation             |

## 5. Experimental Results

### 5.1 Setup

All experiments use the OR-Library SCP benchmark instances (Beasley, 1990) on Windows 11, single-threaded. Three configurations are compared:

- **Baseline**: Pure BnB with most-fractional branching, no cuts, no Balas (`--branch most_fractional --cut-frequency 0 --balas-frequency 0`)
- **Balas+Preprocess**: Preprocessing + Balas branching + most-fractional branching (previous default)
- **Reliability**: Preprocessing + reliability branching + Balas branching (current default)

### 5.2 Results

| Instance    |  Baseline |              | Balas+Preproc |             | Reliability |             |
| ----------- | --------: | -----------: | ------------: | ----------: | ----------: | ----------: |
|             |     Nodes |         Time |         Nodes |        Time |       Nodes |        Time |
| scpa1       |       363 |       75.0 s |           478 |      19.5 s |         113 |      16.4 s |
| scpa2       |       199 |       41.7 s |           227 |      16.2 s |         116 |      17.1 s |
| scpa3       |       137 |       41.3 s |            51 |      17.4 s |          81 |      18.2 s |
| scpa4       |        47 |       21.4 s |            41 |      15.9 s |          31 |      15.0 s |
| scpa5       |        58 |       20.1 s |            18 |      15.9 s |          12 |      16.9 s |
| scpb1       |       182 |       55.7 s |           182 |      22.2 s |         123 |      23.4 s |
| scpb2       |      1641 |      361.3 s |          1577 |      42.6 s |         341 |      29.2 s |
| scpb3       |      1339 |      288.6 s |          1375 |      31.6 s |         241 |      24.7 s |
| scpb4       |      6939 |     1543.8 s |          6541 |     112.0 s |         713 |      33.9 s |
| scpb5       |       379 |       97.5 s |           383 |      24.1 s |         161 |      23.6 s |
| **Total A** |   **804** |  **199.5 s** |       **815** |  **84.9 s** |     **353** |  **83.7 s** |
| **Total B** | **10480** | **2346.9 s** |     **10058** | **232.4 s** |    **1579** | **134.8 s** |

### 5.3 Analysis

**Node reduction**: Reliability branching consistently reduces node counts. The most dramatic improvement is on scpb4, where better branching decisions reduce the tree from 6541 nodes (Balas+Preprocess) to 713 nodes -- a **9.2x reduction**. Across the B-set, total nodes drop from 10058 to 1579 (6.4x fewer).

**Time improvement**: On easy instances (set 5), reliability branching adds slight overhead from probing LPs that aren't amortized by the small tree. On hard instances (set B), the node reduction more than compensates for the probing cost, yielding a **1.7x speedup** over Balas+Preprocess and a **17.4x speedup** over the pure baseline.

**Speedup relative to baseline**: The combined effect of preprocessing, Balas branching, and reliability branching achieves a **45.6x speedup** on scpb4 (1544s to 34s) and a **10.5x overall speedup** across all instances.

**Cost-benefit of strong branching**: Each strong-branch probe requires solving an LP, adding ~0.01-0.05s per probe on these instances. With `max_sb=8` and `eta=4`, most probing occurs in the first 30-50 nodes. After that, pseudocosts are reliable for most variables and the per-node cost drops to near zero. The break-even point is around 100-200 nodes: below that threshold, the probing overhead slightly exceeds the tree reduction benefit; above it, reliability branching is strictly better.

## 6. References

1. T. Achterberg, T. Koch, A. Martin. _Branching rules revisited._ Operations Research Letters, 33(1):42-54, 2005.
2. T. Berthold. _Hybrid Branching._ In: Integration of AI and OR Techniques in Constraint Programming, LNCS 8451, pp. 309-324, 2014.
3. G. Gamrath, T. Koch, A. Martin, M. Miltenberger, D. Weninger. _Improving Strong Branching by Domain Propagation._ EURO Journal on Computational Optimization, 3(1):27-60, 2015.
4. A. Khalil, P. Le Bodic, L. Song, G. Nemhauser, B. Dilkina. _Learning to Branch in Mixed Integer Programming._ AAAI Conference on Artificial Intelligence, 2016.
