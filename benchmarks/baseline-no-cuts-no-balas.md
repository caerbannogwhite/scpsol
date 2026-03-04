# Baseline Benchmark: scpsol --no-cuts --no-balas

Pure branch-and-bound with LP relaxation (HiGHS simplex), preprocessing enabled, no cutting planes, no Balas branching.

**Solver:** scpsol (HiGHS v1.9.0)
**Command:** `scpsol <instance> --no-cuts --no-balas --verbosity 0`
**Date:** 2026-03-04
**Machine:** Windows 11, single-threaded

## Results

| Instance | Primal | Dual | MIP gap | Nodes | LP solves | Wall time |
|----------|-------:|-----:|--------:|------:|----------:|----------:|
| scp51    |    253 |  253 |  0.000% |    12 |         — |    1.08 s |
| scp52    |    302 |  302 |  0.000% |    17 |         — |    1.37 s |
| scp53    |    226 |  226 |  0.000% |     0 |         — |    1.06 s |
| scp54    |    242 |  242 |  0.000% |    15 |         — |    1.16 s |
| scp55    |    211 |  211 |  0.000% |     0 |         — |    1.03 s |
| scp56    |    213 |  213 |  0.000% |     0 |         — |    1.06 s |
| scp57    |    293 |  293 |  0.000% |    27 |         — |    1.21 s |
| scp58    |    288 |  288 |  0.000% |    11 |         — |    1.21 s |
| scp59    |    279 |  279 |  0.000% |     0 |         — |    0.91 s |
| scp510   |    265 |  265 |  0.000% |     0 |         — |    1.05 s |
| scpa1    |    253 |  253 |  0.000% |   363 |         — |   74.98 s |s
| scpa2    |    252 |  252 |  0.000% |   199 |         — |   41.70 s |
| scpa3    |    232 |  232 |  0.000% |   137 |         — |   41.26 s |
| scpa4    |    234 |  234 |  0.000% |    47 |         — |   21.40 s |
| scpa5    |    236 |  236 |  0.000% |    58 |         — |   20.13 s |
| scpb1    |     69 |   69 |  0.000% |   182 |         — |   55.74 s |
| scpb2    |     76 |   76 |  0.000% |  1641 |         — |  361.26 s |
| scpb3    |     80 |   80 |  0.000% |  1339 |         — |  288.56 s |
| scpb4    |     79 |   79 |  0.000% |  6939 |         — | 1543.77 s |
| scpb5    |     72 |   72 |  0.000% |   379 |         — |   97.53 s |

## Summary

- **Set 5** (200×2000): All optimal in ~1 s, most solved at root node.
- **Set A** (300×3000): 20–75 s, up to 363 nodes.
- **Set B** (300×3000, harder costs): 56–1544 s, up to 6939 nodes. scpb4 is the hardest at ~26 min.
- All 20 instances solved to proven optimality (0% MIP gap).
