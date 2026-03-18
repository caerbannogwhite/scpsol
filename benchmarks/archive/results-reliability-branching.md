# Benchmark: Reliability Branching + Preprocessing

Branch-and-bound with LP relaxation (HiGHS simplex), OR-Tools-style preprocessing (essential columns, row domination, probing, reduced-cost fixing), reliability branching with propagation-enhanced strong branching (eta=4, max 8 probes/node), Balas branching (frequency 0.6), no mid-BnB cuts.

**Solver:** scpsol (HiGHS v1.9.0)
**Command:** `scpsol <instance> --verbosity 0` (default settings)
**Date:** 2026-03-11
**Machine:** Windows 11, single-threaded

## Results

| Instance | Primal | Dual | MIP gap | Nodes | Baseline nodes | Wall time | Baseline time | Speedup |
|----------|-------:|-----:|--------:|------:|---------------:|----------:|--------------:|--------:|
| scp41    |    429 |  429 |  0.000% |     0 |            — |    0.24 s |           — |     — |
| scp42    |    512 |  512 |  0.000% |     0 |            — |    0.21 s |           — |     — |
| scp43    |    516 |  516 |  0.000% |     0 |            — |    0.24 s |           — |     — |
| scp44    |    494 |  494 |  0.000% |     2 |            — |    0.16 s |           — |     — |
| scp45    |    512 |  512 |  0.000% |     0 |            — |    0.22 s |           — |     — |
| scp46    |    560 |  560 |  0.000% |     7 |            — |    0.19 s |           — |     — |
| scp47    |    430 |  430 |  0.000% |     0 |            — |    0.19 s |           — |     — |
| scp48    |    492 |  492 |  0.000% |     9 |            — |    0.19 s |           — |     — |
| scp49    |    641 |  641 |  0.000% |    13 |            — |    0.30 s |           — |     — |
| scp410   |    514 |  514 |  0.000% |     3 |            — |    0.14 s |           — |     — |
| scp51    |    253 |  253 |  0.000% |    13 |           12 |    0.94 s |      1.08 s |  1.1x |
| scp52    |    302 |  302 |  0.000% |    42 |           17 |    1.27 s |      1.37 s |  1.1x |
| scp53    |    226 |  226 |  0.000% |     0 |            0 |    1.07 s |      1.06 s |  1.0x |
| scp54    |    242 |  242 |  0.000% |     6 |           15 |    0.91 s |      1.16 s |  1.3x |
| scp55    |    211 |  211 |  0.000% |     0 |            0 |    1.04 s |      1.03 s |  1.0x |
| scp56    |    213 |  213 |  0.000% |     0 |            0 |    0.81 s |      1.06 s |  1.3x |
| scp57    |    293 |  293 |  0.000% |     3 |           27 |    0.86 s |      1.21 s |  1.4x |
| scp58    |    288 |  288 |  0.000% |    10 |           11 |    0.88 s |      1.21 s |  1.4x |
| scp59    |    279 |  279 |  0.000% |     0 |            0 |    0.90 s |      0.91 s |  1.0x |
| scp510   |    265 |  265 |  0.000% |     0 |            0 |    1.04 s |      1.05 s |  1.0x |
| scpa1    |    253 |  253 |  0.000% |   113 |          363 |   16.41 s |     74.98 s |  4.6x |
| scpa2    |    252 |  252 |  0.000% |   116 |          199 |   17.12 s |     41.70 s |  2.4x |
| scpa3    |    232 |  232 |  0.000% |    81 |          137 |   18.24 s |     41.26 s |  2.3x |
| scpa4    |    234 |  234 |  0.000% |    31 |           47 |   15.02 s |     21.40 s |  1.4x |
| scpa5    |    236 |  236 |  0.000% |    12 |           58 |   16.86 s |     20.13 s |  1.2x |
| scpb1    |     69 |   69 |  0.000% |   123 |          182 |   23.41 s |     55.74 s |  2.4x |
| scpb2    |     76 |   76 |  0.000% |   341 |         1641 |   29.19 s |    361.26 s | 12.4x |
| scpb3    |     80 |   80 |  0.000% |   241 |         1339 |   24.68 s |    288.56 s | 11.7x |
| scpb4    |     79 |   79 |  0.000% |   713 |         6939 |   33.86 s |   1543.77 s | 45.6x |
| scpb5    |     72 |   72 |  0.000% |   161 |          379 |   23.63 s |     97.53 s |  4.1x |

## Summary

- **Set 4** (200x1000): All optimal in < 0.3 s (no baseline comparison available).
- **Set 5** (200x2000): All optimal in ~1 s, 1.0-1.4x speedup over baseline.
- **Set A** (300x3000): 15-18 s, 1.2-4.6x speedup, 1.4-3.2x fewer nodes.
- **Set B** (300x3000, harder costs): 23-34 s, 2.4-45.6x speedup, 1.5-9.7x fewer nodes.
- All 30 instances solved to proven optimality (0% MIP gap).
- Total wall time for sets 5+A+B: **243 s** vs 2554 s baseline (**10.5x overall speedup**).
- scpb4 improved from 1544 s (26 min) to 34 s, a **45.6x speedup** with 9.7x fewer nodes.
