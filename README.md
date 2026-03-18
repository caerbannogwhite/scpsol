# scpsol

A Set Covering Problem (SCP) solver built from scratch in C++17, using [HiGHS](https://github.com/ERGO-Code/HiGHS) as the LP backend. Solves 30 standard OR-Library instances (sets 4, 5, A, B) to proven optimality in ~25 seconds total.

## Features

**Preprocessing**
- Iterated fixpoint reduction: essential column fixing, row domination, single/two-column dominance, cost-driven replacement
- Greedy multi-column dominance finder with pipeline reordering
- Budget pruning using incumbent bound
- Bitset-accelerated subset checks (O(m/64) word operations)
- Probing-based row reduction

**Branch-and-Bound**
- Best-bound node selection (priority queue)
- Reliability branching with pseudocost tracking and strong branching with propagation
- Balas additive multi-way branching (stagnation-triggered)
- Per-node LP basis warm-starting, essential column propagation, and reduced-cost fixing
- Mid-BnB cut generation and column removal

**Cuts**
- Chvatal-Gomory cuts (dual-aggregated and row-pair separators)
- Balas cover inequalities (additive procedure)
- Root cut rounds with bound-improvement validation

**Heuristics**
- Chvatal greedy heuristic with redundancy removal (initial incumbent)
- Nearest-integer fixing
- Dual-guided cover repair with incremental coverage

## Benchmark Results

All instances solved to proven optimality with 0% MIP gap (300s time limit, single-threaded):

| Category | Instances | Total Time | Avg Nodes |
|----------|-----------|------------|-----------|
| scp4 (10) | 10/10 | 0.4s | 3 |
| scp5 (10) | 10/10 | 0.6s | 4 |
| scpa (5) | 5/5 | 5.4s | 50 |
| scpb (5) | 5/5 | 19.0s | 262 |
| **Core 30** | **30/30** | **25.4s** | — |

On harder NRE/NRF instances (500 rows, 5000 cols): 7/10 optimal, max residual gap 7.1%.

See [benchmarks/results/comparison.md](benchmarks/results/comparison.md) for the full cross-version comparison (8 tagged versions, 53x cumulative speedup on core instances).

## Build

### Linux / macOS

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

### Windows (Visual Studio)

```bash
cmake -B build
cmake --build build --config Release
```

The binary is placed in `build/Release/scpsol.exe` (Windows) or `build/scpsol` (Linux/macOS).

HiGHS v1.9.0 is fetched automatically via CMake FetchContent.

## Usage

```bash
./scpsol <scp_file> [options]
```

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `--verbosity N` | Verbosity level (0=silent, 1=summary, 2=progress) | `2` |
| `--max-nodes N` | Maximum BnB nodes | `100000` |
| `--time-limit S` | Time limit in seconds | `0` (unlimited) |
| `--branch STRATEGY` | Branch strategy: `reliability`, `most_fractional`, `highest_cost` | `reliability` |
| `--reliability-eta N` | Min observations before trusting pseudocosts | `4` |
| `--reliability-sb N` | Max strong-branch probes per node | `8` |
| `--cut-frequency F` | Mid-BnB cut frequency 0.0-1.0 (0=disable) | `0.0` |
| `--balas-frequency F` | Balas branching frequency 0.0-1.0 (0=disable) | `0.6` |
| `--cut-rounds N` | Mid-BnB cut rounds per event | `3` |
| `--preprocess RULES` | Preprocess rules | `"single,two"` |
| `--show-solution` | Print selected columns | |

### Example

```bash
./scpsol data/scp41.txt --time-limit 60
```

```
Status:    Optimal
Primal:    429
Dual:      429
MIP gap:   0.000000%
Nodes:     0
LP solves: 3
Wall time: 0.014 s
```

## Input Format

Standard SCP text format (OR-Library):

```
nrows ncols
cost_1 cost_2 ... cost_ncols
num_cols_covering_row_1  col_idx_1 col_idx_2 ...
num_cols_covering_row_2  col_idx_1 col_idx_2 ...
...
```

Column indices are 1-based. All constraint coefficients are binary (0/1).

## Test Instances

The `data/` directory contains standard SCP benchmark instances from [OR-Library](http://people.brunel.ac.uk/~mastjjb/jeb/orlib/scpinfo.html) (Beasley) including sets 4, 5, A, B, CLR, CYC, NRE, NRF, NRG, and NRH.

## Docker

```bash
docker build -t scpsol .
docker run --rm scpsol data/scp41.txt
```
