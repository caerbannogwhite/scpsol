# scpsol

A Set Covering Problem (SCP) solver built from scratch in C++17, using [HiGHS](https://github.com/ERGO-Code/HiGHS) as the LP backend.

## Features

- Branch-and-bound with LP relaxation (HiGHS simplex)
- Balas additive branching
- Cutting planes (lifted cover inequalities)
- Primal heuristics (greedy, rounding)
- Preprocessing (singleton/doubleton row reduction)

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

HiGHS v1.9.0 is fetched automatically via CMake FetchContent.

## Usage

```bash
./scpsol <scp_file> [options]
```

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `--verbosity N` | Verbosity level | 2 |
| `--max-nodes N` | Maximum BnB nodes | 100000 |
| `--time-limit S` | Time limit in seconds | 0 (unlimited) |
| `--no-cuts` | Disable cutting planes | |
| `--no-balas` | Disable Balas branching | |
| `--show-solution` | Print selected columns | |
| `--preprocess RULES` | Preprocess rules | `"single,two"` |

### Example

```bash
./scpsol data/scp41.txt --time-limit 60
```

```
Status:    Optimal
Primal:    429
Dual:      429
MIP gap:   0.000000%
Nodes:     34
LP solves: 87
Wall time: 0.412 s
```

## Input format

Standard SCP text format (OR-Library):

```
nrows ncols
cost_1 cost_2 ... cost_ncols
num_cols_covering_row_1  col_idx_1 col_idx_2 ...
num_cols_covering_row_2  col_idx_1 col_idx_2 ...
...
```

Column indices are 1-based. All constraint coefficients are binary (0/1).

## Test instances

The `data/` directory contains standard SCP benchmark instances from OR-Library (Beasley) including sets 4, 5, A, B, CLR, CYC, NRE, NRF, NRG, and NRH.

## Docker

```bash
docker build -t scpsol .
docker run --rm scpsol data/scp41.txt
```
