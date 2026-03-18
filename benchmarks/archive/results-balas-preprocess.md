# Benchmark: Balas + Preprocessing

Branch-and-bound with LP relaxation (HiGHS simplex), OR-Tools-style preprocessing (essential columns, row domination, probing, reduced-cost fixing), Balas branching (frequency 0.6), no mid-BnB cuts.

**Solver:** scpsol (HiGHS v1.9.0)
**Command:** `scpsol <instance> --verbosity 0` (default settings: `--cut-frequency 0 --balas-frequency 0.6`)
**Date:** 2026-03-11
**Machine:** Windows 11, single-threaded

## Results

| Instance | Primal | Dual | MIP gap | Nodes | Wall time | Baseline time | Speedup |
|----------|-------:|-----:|--------:|------:|----------:|--------------:|--------:|
| scp41    |    429 |  429 |  0.000% |     0 |    0.24 s |             — |       — |
| scp42    |    512 |  512 |  0.000% |     0 |    0.23 s |             — |       — |
| scp43    |    516 |  516 |  0.000% |     0 |    0.24 s |             — |       — |
| scp44    |    494 |  494 |  0.000% |     6 |    0.17 s |             — |       — |
| scp45    |    512 |  512 |  0.000% |     0 |    0.24 s |             — |       — |
| scp46    |    560 |  560 |  0.000% |    17 |    0.18 s |             — |       — |
| scp47    |    430 |  430 |  0.000% |     0 |    0.20 s |             — |       — |
| scp48    |    492 |  492 |  0.000% |    22 |    0.17 s |             — |       — |
| scp49    |    641 |  641 |  0.000% |    12 |    0.22 s |             — |       — |
| scp410   |    514 |  514 |  0.000% |     2 |    0.13 s |             — |       — |
| scp51    |    253 |  253 |  0.000% |    52 |    0.90 s |        1.08 s |    1.2x |
| scp52    |    302 |  302 |  0.000% |    35 |    0.94 s |        1.37 s |    1.5x |
| scp53    |    226 |  226 |  0.000% |     0 |    1.06 s |        1.06 s |    1.0x |
| scp54    |    242 |  242 |  0.000% |    37 |    0.88 s |        1.16 s |    1.3x |
| scp55    |    211 |  211 |  0.000% |     0 |    1.03 s |        1.03 s |    1.0x |
| scp56    |    213 |  213 |  0.000% |     0 |    0.81 s |        1.06 s |    1.3x |
| scp57    |    293 |  293 |  0.000% |     8 |    0.84 s |        1.21 s |    1.4x |
| scp58    |    288 |  288 |  0.000% |     4 |    0.81 s |        1.21 s |    1.5x |
| scp59    |    279 |  279 |  0.000% |     0 |    0.91 s |        0.91 s |    1.0x |
| scp510   |    265 |  265 |  0.000% |     0 |    1.06 s |        1.05 s |    1.0x |
| scpa1    |    253 |  253 |  0.000% |   478 |   19.49 s |       74.98 s |    3.8x |
| scpa2    |    252 |  252 |  0.000% |   227 |   16.15 s |       41.70 s |    2.6x |
| scpa3    |    232 |  232 |  0.000% |    51 |   17.37 s |       41.26 s |    2.4x |
| scpa4    |    234 |  234 |  0.000% |    41 |   15.88 s |       21.40 s |    1.3x |
| scpa5    |    236 |  236 |  0.000% |    18 |   15.87 s |       20.13 s |    1.3x |
| scpb1    |     69 |   69 |  0.000% |   182 |   22.22 s |       55.74 s |    2.5x |
| scpb2    |     76 |   76 |  0.000% |  1577 |   42.56 s |      361.26 s |    8.5x |
| scpb3    |     80 |   80 |  0.000% |  1375 |   31.58 s |      288.56 s |    9.1x |
| scpb4    |     79 |   79 |  0.000% |  6541 |  111.98 s |     1543.77 s |   13.8x |
| scpb5    |     72 |   72 |  0.000% |   383 |   24.10 s |       97.53 s |    4.0x |

## Summary

- **Set 4** (200x1000): All optimal in < 0.25 s (no baseline comparison available).
- **Set 5** (200x2000): All optimal in ~1 s, 1.0-1.5x speedup over baseline.
- **Set A** (300x3000): 16-20 s, 1.3-3.8x speedup over baseline.
- **Set B** (300x3000, harder costs): 22-112 s, 2.5-13.8x speedup over baseline.
- All 30 instances solved to proven optimality (0% MIP gap).
- Total wall time for sets 5+A+B: 326 s vs 2554 s baseline (**7.8x overall speedup**).
