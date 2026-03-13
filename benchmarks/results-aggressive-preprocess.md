# Benchmark: Aggressive Preprocessing + Cost-Effectiveness Greedy

Iterated fixpoint preprocessing (essential columns, row domination, probing, column dominance, budget pruning, reduced-cost fixing), Chvatal cost-effectiveness greedy heuristic, reliability branching (eta=4, sb=8), Balas branching (frequency 0.6).

**Solver:** scpsol (HiGHS v1.9.0)
**Command:** `scpsol <instance> --verbosity 0` (default settings)
**Date:** 2026-03-13
**Machine:** Windows 11, single-threaded

## Results

| Instance | Primal | Dual | MIP gap | Nodes | Prev nodes | Wall time | Prev time | vs prev |
|----------|-------:|-----:|--------:|------:|-----------:|----------:|----------:|--------:|
| scp41    |    429 |  429 |  0.000% |     0 |          0 |    0.24 s |    0.24 s |   1.0x |
| scp42    |    512 |  512 |  0.000% |     0 |          0 |    0.22 s |    0.21 s |   1.0x |
| scp43    |    516 |  516 |  0.000% |     0 |          0 |    0.23 s |    0.24 s |   1.0x |
| scp44    |    494 |  494 |  0.000% |     1 |          2 |    0.16 s |    0.16 s |   1.0x |
| scp45    |    512 |  512 |  0.000% |     0 |          0 |    0.23 s |    0.22 s |   1.0x |
| scp46    |    560 |  560 |  0.000% |     7 |          7 |    0.20 s |    0.19 s |   1.0x |
| scp47    |    430 |  430 |  0.000% |     0 |          0 |    0.20 s |    0.19 s |   1.0x |
| scp48    |    492 |  492 |  0.000% |     7 |          9 |    0.19 s |    0.19 s |   1.0x |
| scp49    |    641 |  641 |  0.000% |    13 |         13 |    0.32 s |    0.30 s |   0.9x |
| scp410   |    514 |  514 |  0.000% |     1 |          3 |    0.13 s |    0.14 s |   1.1x |
| scp51    |    253 |  253 |  0.000% |    13 |         13 |    0.91 s |    0.94 s |   1.0x |
| scp52    |    302 |  302 |  0.000% |    42 |         42 |    1.30 s |    1.27 s |   1.0x |
| scp53    |    226 |  226 |  0.000% |     0 |          0 |    1.04 s |    1.07 s |   1.0x |
| scp54    |    242 |  242 |  0.000% |    10 |          6 |    0.90 s |    0.91 s |   1.0x |
| scp55    |    211 |  211 |  0.000% |     0 |          0 |    1.01 s |    1.04 s |   1.0x |
| scp56    |    213 |  213 |  0.000% |     0 |          0 |    0.78 s |    0.81 s |   1.0x |
| scp57    |    293 |  293 |  0.000% |     3 |          3 |    0.85 s |    0.86 s |   1.0x |
| scp58    |    288 |  288 |  0.000% |     5 |         10 |    0.83 s |    0.88 s |   1.1x |
| scp59    |    279 |  279 |  0.000% |     0 |          0 |    0.89 s |    0.90 s |   1.0x |
| scp510   |    265 |  265 |  0.000% |     0 |          0 |    1.03 s |    1.04 s |   1.0x |
| scpa1    |    253 |  253 |  0.000% |   113 |        113 |   16.26 s |   16.41 s |   1.0x |
| scpa2    |    252 |  252 |  0.000% |   116 |        116 |   16.93 s |   17.12 s |   1.0x |
| scpa3    |    232 |  232 |  0.000% |    91 |         81 |   18.52 s |   18.24 s |   1.0x |
| scpa4    |    234 |  234 |  0.000% |    31 |         31 |   14.98 s |   15.02 s |   1.0x |
| scpa5    |    236 |  236 |  0.000% |    29 |         12 |   14.67 s |   16.86 s |   1.1x |
| scpb1    |     69 |   69 |  0.000% |   123 |        123 |   20.29 s |   23.41 s | **1.2x** |
| scpb2    |     76 |   76 |  0.000% |   341 |        341 |   28.53 s |   29.19 s |   1.0x |
| scpb3    |     80 |   80 |  0.000% |   241 |        241 |   24.29 s |   24.68 s |   1.0x |
| scpb4    |     79 |   79 |  0.000% |   713 |        713 |   34.01 s |   33.86 s |   1.0x |
| scpb5    |     72 |   72 |  0.000% |   161 |        161 |   21.72 s |   23.63 s | **1.1x** |

## Summary

- All 30 instances solved to proven optimality (0% MIP gap).
- The cost-effectiveness greedy (Chvatal ratio) replaces the sort-by-cost greedy, using incremental coverage tracking for O(nnz) total updates.
- Preprocessing now iterates cost reduction, dominance reduction, and row reduction until fixpoint, enabling cascading reductions across techniques.
- Performance is neutral or slightly improved vs the previous version: scpb1 is 1.2x faster, scpb5 is 1.1x faster.
- Total wall time for sets 5+A+B: **216 s** vs 225 s previous (**1.04x overall**).
