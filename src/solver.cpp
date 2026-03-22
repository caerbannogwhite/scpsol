#include "solver.h"
#include "balas.h"
#include "bnb.h"
#include "cuts.h"
#include "decomposition.h"
#include "diving.h"
#include "heuristics.h"
#include "rins.h"
#include "lagrangian.h"
#include "lp.h"
#include "preprocessor.h"
#include "reliability.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <random>
#include <unordered_map>
#include <memory>
#include <queue>
#include <string>

namespace scpsol {

// ---- Symmetry breaking: orbital fixing ----
// Columns with the same cost that cover the same set of *undecided* rows
// at a BnB node are interchangeable. We keep one representative per
// equivalence class and fix the rest to 0.

struct OrbitalFixingState {
    // Column groups with identical cost (only groups of size > 1)
    std::vector<std::vector<int>> cost_groups;
    // Zobrist hash table: row index -> random 64-bit value
    std::vector<uint64_t> row_hash;
    int total_fixed = 0;
    int calls = 0;
    bool disabled = false; // auto-disable after fruitless attempts
    // Scratch buffers (reused to avoid per-call allocation)
    std::vector<int8_t> col_state_buf;
    std::vector<uint8_t> row_decided_buf;
};

static OrbitalFixingState build_orbital_state(
    const BaseRelaxationModel &base) {
    OrbitalFixingState state;

    // Build Zobrist hash table
    std::mt19937_64 rng(0xBEEF5EED);
    state.row_hash.resize(static_cast<size_t>(base.nrows));
    for (int i = 0; i < base.nrows; ++i)
        state.row_hash[static_cast<size_t>(i)] = rng();

    // Group columns by cost (quantized to avoid floating-point issues)
    std::unordered_map<int64_t, std::vector<int>> cost_map;
    for (int j = 0; j < base.ncols; ++j) {
        int64_t key = static_cast<int64_t>(std::round(base.obj[static_cast<size_t>(j)] * 1e6));
        cost_map[key].push_back(j);
    }
    for (auto &[key, cols] : cost_map) {
        if (cols.size() > 1)
            state.cost_groups.push_back(std::move(cols));
    }
    return state;
}

// Compute orbital fixings at a BnB node.
// Returns additional decisions (fix_value=0) for equivalent columns.
static std::vector<BranchDecision> compute_orbital_fixings(
    const ScpAdjacency &adj,
    const std::vector<BranchDecision> &decisions,
    OrbitalFixingState &state) {

    std::vector<BranchDecision> fixings;
    if (state.cost_groups.empty() || state.disabled) return fixings;

    const int ncols = adj.ncols;
    const int nrows = adj.nrows;

    // Reuse scratch buffers
    state.col_state_buf.assign(static_cast<size_t>(ncols), 0);
    state.row_decided_buf.assign(static_cast<size_t>(nrows), 0);
    auto &col_state = state.col_state_buf;
    auto &row_decided = state.row_decided_buf;

    for (const auto &d : decisions) {
        if (d.var_index < 0 || d.var_index >= ncols) continue;
        col_state[static_cast<size_t>(d.var_index)] =
            static_cast<int8_t>(d.fix_value == 1 ? 1 : -1);
        if (d.fix_value == 1) {
            for (int r : adj.rows_by_col[static_cast<size_t>(d.var_index)])
                if (r >= 0 && r < nrows)
                    row_decided[static_cast<size_t>(r)] = 1;
        }
    }

    // For each cost group, find equivalence classes of free columns
    // based on their residual coverage (undecided rows).
    for (const auto &group : state.cost_groups) {
        // Hash free columns by their residual coverage
        struct HashEntry { uint64_t hash; int col; };
        std::vector<HashEntry> entries;
        entries.reserve(group.size());
        for (int j : group) {
            if (col_state[static_cast<size_t>(j)] != 0) continue;
            uint64_t h = 0;
            int cnt = 0;
            for (int r : adj.rows_by_col[static_cast<size_t>(j)]) {
                if (r >= 0 && r < nrows && !row_decided[static_cast<size_t>(r)]) {
                    h ^= state.row_hash[static_cast<size_t>(r)];
                    ++cnt;
                }
            }
            // Mix count into hash to reduce collisions
            h ^= static_cast<uint64_t>(cnt) * 0x9E3779B97F4A7C15ULL;
            entries.push_back({h, j});
        }
        if (entries.size() <= 1) continue;

        // Sort by hash to cluster potential equivalences
        std::sort(entries.begin(), entries.end(),
                  [](const HashEntry &a, const HashEntry &b) {
                      return a.hash < b.hash;
                  });

        // Process each hash bucket
        size_t i = 0;
        while (i < entries.size()) {
            size_t j = i + 1;
            while (j < entries.size() && entries[j].hash == entries[i].hash)
                ++j;
            // entries[i..j) have the same hash
            if (j - i > 1) {
                // Verify exact equivalence: build residual row sets and compare
                // Use first column as reference
                const int ref_col = entries[i].col;
                std::vector<int> ref_residual;
                for (int r : adj.rows_by_col[static_cast<size_t>(ref_col)])
                    if (r >= 0 && r < nrows && !row_decided[static_cast<size_t>(r)])
                        ref_residual.push_back(r);

                for (size_t k = i + 1; k < j; ++k) {
                    const int cand = entries[k].col;
                    // Build candidate residual and compare
                    bool match = true;
                    const auto &cand_rows = adj.rows_by_col[static_cast<size_t>(cand)];
                    size_t ri = 0, ci = 0;
                    // Both row lists are sorted; compare via merge
                    const auto &ref_rows = adj.rows_by_col[static_cast<size_t>(ref_col)];
                    ri = 0; ci = 0;
                    size_t rr = 0; // index into ref_residual for validation
                    // Rebuild residual for candidate inline and compare
                    std::vector<int> cand_residual;
                    for (int r : cand_rows)
                        if (r >= 0 && r < nrows && !row_decided[static_cast<size_t>(r)])
                            cand_residual.push_back(r);
                    if (cand_residual != ref_residual)
                        match = false;

                    if (match) {
                        fixings.push_back({cand, 0});
                        col_state[static_cast<size_t>(cand)] = -1;
                    }
                }
            }
            i = j;
        }
    }
    return fixings;
}

static void adopt_incumbent_solution(
    std::vector<double> &best_solution,
    const std::vector<double> &active_solution,
    int ncols,
    int ncols_input,
    const std::vector<int> &active_to_original) {
    best_solution.assign(static_cast<size_t>(ncols_input), 0.0);
    const int copy_cols = std::min(ncols, static_cast<int>(active_solution.size()));
    for (int j = 0; j < copy_cols; ++j) {
        if (active_solution[static_cast<size_t>(j)] > 0.5) {
            const int orig = active_to_original[static_cast<size_t>(j)];
            if (orig >= 0 && orig < ncols_input) {
                best_solution[static_cast<size_t>(orig)] = 1.0;
            }
        }
    }
}

// Best-bound frontier: nodes ordered by ascending parent_dual_bound.
// Uses a vector managed as a min-heap for efficient best-bound selection.
struct BestBoundFrontier {
    std::vector<int> heap;
    const std::vector<BranchNodeState> *nodes_ptr = nullptr;
    double min_raw_bound = std::numeric_limits<double>::infinity();

    void init(const std::vector<BranchNodeState> *np) { nodes_ptr = np; }

    bool empty() const { return heap.empty(); }
    size_t size() const { return heap.size(); }

    // Tightened dual bound of the best node (heap minimum)
    double top_bound() const {
        return (*nodes_ptr)[static_cast<size_t>(heap[0])].parent_dual_bound;
    }

    bool cmp(int a, int b) const {
        return (*nodes_ptr)[static_cast<size_t>(a)].parent_dual_bound >
               (*nodes_ptr)[static_cast<size_t>(b)].parent_dual_bound;
    }

    void push(int idx) {
        heap.push_back(idx);
        std::push_heap(heap.begin(), heap.end(),
            [this](int a, int b) { return cmp(a, b); });
        const double raw = (*nodes_ptr)[static_cast<size_t>(idx)].parent_dual_bound_raw;
        if (raw < min_raw_bound) min_raw_bound = raw;
    }

    int pop() {
        std::pop_heap(heap.begin(), heap.end(),
            [this](int a, int b) { return cmp(a, b); });
        int idx = heap.back();
        heap.pop_back();
        // If the popped node held the minimum raw bound, recompute
        const double popped_raw = (*nodes_ptr)[static_cast<size_t>(idx)].parent_dual_bound_raw;
        if (popped_raw <= min_raw_bound + 1e-12)
            recompute_min_raw();
        return idx;
    }

    void prune(double best_obj, double tol, int verbosity,
               std::unordered_map<int, HighsBasis> *bases = nullptr) {
        const double prune_bound = best_obj - tol;
        size_t before = heap.size();
        std::vector<int> surviving;
        surviving.reserve(heap.size());
        min_raw_bound = std::numeric_limits<double>::infinity();
        for (const int idx : heap) {
            const auto &nd = (*nodes_ptr)[static_cast<size_t>(idx)];
            if (nd.parent_dual_bound < prune_bound) {
                surviving.push_back(idx);
                if (nd.parent_dual_bound_raw < min_raw_bound)
                    min_raw_bound = nd.parent_dual_bound_raw;
            } else if (bases) {
                bases->erase(idx);
            }
        }
        heap = std::move(surviving);
        std::make_heap(heap.begin(), heap.end(),
            [this](int a, int b) { return cmp(a, b); });
        if (heap.size() < before && verbosity >= 3) {
            fprintf(stderr, "           Frontier pruned: %zu -> %zu nodes\n", before, heap.size());
        }
    }

    void rebuild_heap() {
        std::make_heap(heap.begin(), heap.end(),
            [this](int a, int b) { return cmp(a, b); });
    }

    void recompute_min_raw() {
        min_raw_bound = std::numeric_limits<double>::infinity();
        for (const int idx : heap) {
            const double raw = (*nodes_ptr)[static_cast<size_t>(idx)].parent_dual_bound_raw;
            if (raw < min_raw_bound) min_raw_bound = raw;
        }
    }
};

static int mid_bnb_column_removal(
    BaseRelaxationModel &base,
    double best_obj,
    double tol,
    BestBoundFrontier &frontier,
    std::vector<BranchNodeState> &nodes,
    int verbosity,
    PseudocostState *pc_state = nullptr) {
    ModelReductionResult reduction = reduce_base_model(base, best_obj, tol);
    if (reduction.columns_removed <= 0) return 0;
    base.build_transpose();
    if (verbosity >= 3)
        fprintf(stderr, "           Mid-BnB reduction: %d cols removed, %d remaining\n",
                reduction.columns_removed, base.ncols);
    std::vector<int> surviving;
    for (const int idx : frontier.heap) {
        if (remap_branch_node(nodes[static_cast<size_t>(idx)], reduction.old_to_new))
            surviving.push_back(idx);
    }
    frontier.heap = std::move(surviving);
    frontier.rebuild_heap();
    frontier.recompute_min_raw();
    if (pc_state) pc_state->remap(reduction.old_to_new, base.ncols);
    return reduction.columns_removed;
}

static int mid_bnb_budget_pruning(
    BaseRelaxationModel &base,
    double best_obj,
    double tol,
    double preprocess_time_limit,
    BestBoundFrontier &frontier,
    std::vector<BranchNodeState> &nodes,
    int verbosity,
    PseudocostState *pc_state = nullptr) {
    ModelReductionResult reduction = reduce_base_model_budget_pruning(
        base, best_obj, tol, preprocess_time_limit);
    if (reduction.columns_removed <= 0) return 0;
    base.build_transpose();
    if (verbosity >= 3)
        fprintf(stderr, "           Mid-BnB budget pruning: %d cols removed, %d remaining\n",
                reduction.columns_removed, base.ncols);
    std::vector<int> surviving;
    for (const int idx : frontier.heap) {
        if (remap_branch_node(nodes[static_cast<size_t>(idx)], reduction.old_to_new))
            surviving.push_back(idx);
    }
    frontier.heap = std::move(surviving);
    frontier.rebuild_heap();
    frontier.recompute_min_raw();
    if (pc_state) pc_state->remap(reduction.old_to_new, base.ncols);
    return reduction.columns_removed;
}

static void apply_dominance_reduction(
    int nrows, int &ncols,
    std::vector<int> &csr_inds, std::vector<int> &csr_offs, std::vector<double> &csr_vals,
    std::vector<double> &obj,
    std::vector<int> &active_to_input,
    const std::string &rules_config,
    double tol, double time_limit_sec, int verbosity) {
    if (ncols <= 0) return;

    ColumnPreprocessContext ctx;
    ctx.nrows = nrows;
    ctx.ncols = ncols;
    ctx.costs.assign(obj.begin(), obj.begin() + ncols);
    ctx.active.assign(static_cast<size_t>(ncols), 1);
    if (time_limit_sec > 0.0) {
        ctx.deadline = std::chrono::steady_clock::now() +
                       std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                           std::chrono::duration<double>(time_limit_sec));
    }

    ctx.rows_by_column.assign(static_cast<size_t>(ncols), std::vector<int>());
    for (int i = 0; i < nrows; ++i) {
        for (int k = csr_offs[static_cast<size_t>(i)]; k < csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int col = csr_inds[static_cast<size_t>(k)];
            if (col >= 0 && col < ncols && csr_vals[static_cast<size_t>(k)] > tol)
                ctx.rows_by_column[static_cast<size_t>(col)].push_back(i);
        }
    }

    auto rules = make_preprocess_rules(rules_config);
    int total_removed = 0;
    for (const auto &rule : rules) {
        total_removed += rule->apply(ctx, tol);
    }
    if (total_removed <= 0) return;

    // Rebuild
    std::vector<int> old_to_new(static_cast<size_t>(ncols));
    std::vector<int> new_to_old;
    std::vector<int> new_active;
    int new_col = 0;
    for (int j = 0; j < ncols; ++j) {
        if (!ctx.active[static_cast<size_t>(j)]) {
            old_to_new[static_cast<size_t>(j)] = -1;
        } else {
            old_to_new[static_cast<size_t>(j)] = new_col;
            new_to_old.push_back(j);
            new_active.push_back(active_to_input[static_cast<size_t>(j)]);
            ++new_col;
        }
    }

    std::vector<double> new_obj(static_cast<size_t>(new_col));
    for (int j = 0; j < new_col; ++j)
        new_obj[static_cast<size_t>(j)] = obj[static_cast<size_t>(new_to_old[static_cast<size_t>(j)])];

    std::vector<int> new_inds;
    std::vector<int> new_offs;
    std::vector<double> new_vals;
    new_offs.push_back(0);
    for (int i = 0; i < nrows; ++i) {
        for (int k = csr_offs[static_cast<size_t>(i)]; k < csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int c = csr_inds[static_cast<size_t>(k)];
            if (c >= 0 && c < ncols) {
                const int m = old_to_new[static_cast<size_t>(c)];
                if (m >= 0) {
                    new_inds.push_back(m);
                    new_vals.push_back(csr_vals[static_cast<size_t>(k)]);
                }
            }
        }
        new_offs.push_back(static_cast<int>(new_vals.size()));
    }

    if (verbosity >= 3)
        fprintf(stderr, "  Dominance reduction: %d -> %d cols\n", ncols, new_col);

    ncols = new_col;
    csr_inds = std::move(new_inds);
    csr_offs = std::move(new_offs);
    csr_vals = std::move(new_vals);
    obj = std::move(new_obj);
    active_to_input = std::move(new_active);
}

static void apply_cost_reduction(
    int nrows, int &ncols,
    std::vector<int> &csr_inds, std::vector<int> &csr_offs, std::vector<double> &csr_vals,
    std::vector<double> &obj,
    std::vector<int> &active_to_input,
    double incumbent_bound, double tol, int verbosity) {
    if (!std::isfinite(incumbent_bound) || ncols <= 0) return;

    std::vector<int> old_to_new(static_cast<size_t>(ncols));
    std::vector<int> new_to_old;
    std::vector<int> new_active;
    int new_col = 0;
    int removed = 0;
    for (int j = 0; j < ncols; ++j) {
        if (obj[static_cast<size_t>(j)] + tol >= incumbent_bound) {
            old_to_new[static_cast<size_t>(j)] = -1;
            ++removed;
        } else {
            old_to_new[static_cast<size_t>(j)] = new_col;
            new_to_old.push_back(j);
            new_active.push_back(active_to_input[static_cast<size_t>(j)]);
            ++new_col;
        }
    }
    if (removed == 0) return;

    std::vector<double> new_obj(static_cast<size_t>(new_col));
    for (int j = 0; j < new_col; ++j)
        new_obj[static_cast<size_t>(j)] = obj[static_cast<size_t>(new_to_old[static_cast<size_t>(j)])];

    std::vector<int> new_inds;
    std::vector<int> new_offs;
    std::vector<double> new_vals;
    new_offs.push_back(0);
    for (int i = 0; i < nrows; ++i) {
        for (int k = csr_offs[static_cast<size_t>(i)]; k < csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int c = csr_inds[static_cast<size_t>(k)];
            if (c >= 0 && c < ncols) {
                const int m = old_to_new[static_cast<size_t>(c)];
                if (m >= 0) {
                    new_inds.push_back(m);
                    new_vals.push_back(csr_vals[static_cast<size_t>(k)]);
                }
            }
        }
        new_offs.push_back(static_cast<int>(new_vals.size()));
    }

    if (verbosity >= 3)
        fprintf(stderr, "  Cost reduction: %d cols removed, %d remaining\n", removed, new_col);

    ncols = new_col;
    csr_inds = std::move(new_inds);
    csr_offs = std::move(new_offs);
    csr_vals = std::move(new_vals);
    obj = std::move(new_obj);
    active_to_input = std::move(new_active);
}

static void apply_budget_pruning_preprocess(
    int nrows, int &ncols,
    std::vector<int> &csr_inds, std::vector<int> &csr_offs, std::vector<double> &csr_vals,
    std::vector<double> &obj,
    std::vector<int> &active_to_input,
    double incumbent_bound, double tol, double time_limit, int verbosity) {
    if (!std::isfinite(incumbent_bound) || ncols <= 0) return;

    ColumnPreprocessContext ctx;
    ctx.nrows = nrows;
    ctx.ncols = ncols;
    ctx.costs.assign(obj.begin(), obj.begin() + ncols);
    ctx.active.assign(static_cast<size_t>(ncols), 1);
    ctx.incumbent_bound = incumbent_bound;
    if (time_limit > 0.0) {
        ctx.deadline = std::chrono::steady_clock::now() +
                       std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                           std::chrono::duration<double>(time_limit));
    }

    ctx.rows_by_column.assign(static_cast<size_t>(ncols), std::vector<int>());
    for (int i = 0; i < nrows; ++i) {
        for (int k = csr_offs[static_cast<size_t>(i)]; k < csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int col = csr_inds[static_cast<size_t>(k)];
            if (col >= 0 && col < ncols && csr_vals[static_cast<size_t>(k)] > tol)
                ctx.rows_by_column[static_cast<size_t>(col)].push_back(i);
        }
    }

    auto rules = make_preprocess_rules("incumbent_budget");
    int total_removed = 0;
    for (const auto &rule : rules) total_removed += rule->apply(ctx, tol);
    if (total_removed <= 0) return;

    std::vector<int> old_to_new(static_cast<size_t>(ncols));
    std::vector<int> new_to_old;
    std::vector<int> new_active;
    int new_col = 0;
    for (int j = 0; j < ncols; ++j) {
        if (!ctx.active[static_cast<size_t>(j)]) {
            old_to_new[static_cast<size_t>(j)] = -1;
        } else {
            old_to_new[static_cast<size_t>(j)] = new_col;
            new_to_old.push_back(j);
            new_active.push_back(active_to_input[static_cast<size_t>(j)]);
            ++new_col;
        }
    }

    std::vector<double> new_obj(static_cast<size_t>(new_col));
    for (int j = 0; j < new_col; ++j)
        new_obj[static_cast<size_t>(j)] = obj[static_cast<size_t>(new_to_old[static_cast<size_t>(j)])];

    std::vector<int> new_inds, new_offs;
    std::vector<double> new_vals;
    new_offs.push_back(0);
    for (int i = 0; i < nrows; ++i) {
        for (int k = csr_offs[static_cast<size_t>(i)]; k < csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int c = csr_inds[static_cast<size_t>(k)];
            if (c >= 0 && c < ncols) {
                const int m = old_to_new[static_cast<size_t>(c)];
                if (m >= 0) {
                    new_inds.push_back(m);
                    new_vals.push_back(csr_vals[static_cast<size_t>(k)]);
                }
            }
        }
        new_offs.push_back(static_cast<int>(new_vals.size()));
    }

    if (verbosity >= 3)
        fprintf(stderr, "  Budget pruning: %d cols removed, %d remaining\n", total_removed, new_col);

    ncols = new_col;
    csr_inds = std::move(new_inds);
    csr_offs = std::move(new_offs);
    csr_vals = std::move(new_vals);
    obj = std::move(new_obj);
    active_to_input = std::move(new_active);
}

SolverResult solve(const ScpInstance &instance, const SolverConfig &config) {
    using Clock = std::chrono::steady_clock;
    const auto start_time = Clock::now();

    SolverResult result;
    const int verbosity = config.verbosity;
    const double tol = config.feasibility_tol;
    const double integ_tol = config.integrality_tol;

    // Working copies for preprocessing
    int nrows = instance.nrows;
    int ncols = instance.ncols;
    const int ncols_input = ncols;
    std::vector<int> csr_inds = instance.csr_indices;
    std::vector<int> csr_offs = instance.csr_offsets;
    std::vector<double> csr_vals = instance.csr_values;
    std::vector<double> obj = instance.costs;
    std::vector<int> active_to_input(static_cast<size_t>(ncols));
    for (int j = 0; j < ncols; ++j) active_to_input[static_cast<size_t>(j)] = j;

    const bool obj_is_integral = has_integer_objective(obj, ncols, integ_tol);
    if (obj_is_integral && verbosity >= 2)
        fprintf(stderr, "Objective coefficients are integral; enabling dual bound tightening\n");

    double best_obj = std::numeric_limits<double>::infinity();
    std::vector<double> best_solution;
    std::string incumbent_source = "none";
    double global_dual_bound = -std::numeric_limits<double>::infinity();
    double global_dual_bound_raw = -std::numeric_limits<double>::infinity();
    double fixed_preprocess_cost = 0.0;
    std::vector<int> fixed_original_cols;

    auto heuristics = make_integer_heuristics(config.heuristic_config);

    // ================================================================
    // Phase 1: Greedy heuristic
    // ================================================================
    if (verbosity >= 2) fprintf(stderr, "Phase 1: Greedy set cover heuristic\n");
    {
        auto greedy = greedy_set_cover_heuristic(nrows, ncols, csr_inds, csr_offs, csr_vals, obj.data());
        if (greedy.feasible) {
            best_obj = greedy.objective;
            best_solution.assign(static_cast<size_t>(ncols_input), 0.0);
            for (int col : greedy.selected_columns) {
                int input_col = active_to_input[static_cast<size_t>(col)];
                if (input_col >= 0 && input_col < ncols_input)
                    best_solution[static_cast<size_t>(input_col)] = 1.0;
            }
            incumbent_source = "greedy";
            if (verbosity >= 2)
                fprintf(stderr, "  Greedy incumbent: %.12g\n", best_obj);
        }
    }

    // ================================================================
    // Phase 2: Cost + budget + dominance + row reduction (iterated)
    // ================================================================
    // Helper: run Dominance Finder with CSR compaction
    auto run_dominance_finder = [&]() {
        if (ncols <= 0 || nrows <= 0) return;
        ColumnPreprocessContext df_ctx;
        df_ctx.nrows = nrows;
        df_ctx.ncols = ncols;
        df_ctx.costs.resize(static_cast<size_t>(ncols));
        df_ctx.active.assign(static_cast<size_t>(ncols), 1);
        df_ctx.rows_by_column.resize(static_cast<size_t>(ncols));
        df_ctx.deadline = std::chrono::steady_clock::now() +
            std::chrono::seconds(static_cast<int>(config.preprocess_time_limit));

        for (int j = 0; j < ncols; ++j)
            df_ctx.costs[static_cast<size_t>(j)] = obj[static_cast<size_t>(j)];
        for (int i = 0; i < nrows; ++i) {
            for (int k = csr_offs[static_cast<size_t>(i)]; k < csr_offs[static_cast<size_t>(i) + 1]; ++k) {
                const int col = csr_inds[static_cast<size_t>(k)];
                if (col >= 0 && col < ncols)
                    df_ctx.rows_by_column[static_cast<size_t>(col)].push_back(i);
            }
        }

        int df_removed = dominance_finder(df_ctx, tol, 0.5, verbosity);
        if (df_removed > 0) {
            std::vector<int> old_to_new(static_cast<size_t>(ncols), -1);
            int new_idx = 0;
            for (int j = 0; j < ncols; ++j) {
                if (df_ctx.active[static_cast<size_t>(j)])
                    old_to_new[static_cast<size_t>(j)] = new_idx++;
            }
            std::vector<double> new_obj;
            std::vector<int> new_active;
            for (int j = 0; j < ncols; ++j) {
                if (df_ctx.active[static_cast<size_t>(j)]) {
                    new_obj.push_back(obj[static_cast<size_t>(j)]);
                    new_active.push_back(active_to_input[static_cast<size_t>(j)]);
                }
            }
            std::vector<int> new_inds;
            std::vector<int> new_offs = {0};
            std::vector<double> new_vals;
            for (int i = 0; i < nrows; ++i) {
                for (int k = csr_offs[static_cast<size_t>(i)]; k < csr_offs[static_cast<size_t>(i) + 1]; ++k) {
                    const int old_col = csr_inds[static_cast<size_t>(k)];
                    if (old_col >= 0 && old_col < ncols && old_to_new[static_cast<size_t>(old_col)] >= 0) {
                        new_inds.push_back(old_to_new[static_cast<size_t>(old_col)]);
                        new_vals.push_back(csr_vals[static_cast<size_t>(k)]);
                    }
                }
                new_offs.push_back(static_cast<int>(new_inds.size()));
            }
            ncols = new_idx;
            obj = std::move(new_obj);
            active_to_input = std::move(new_active);
            csr_inds = std::move(new_inds);
            csr_offs = std::move(new_offs);
            csr_vals = std::move(new_vals);
            if (verbosity >= 2)
                fprintf(stderr, "  Dominance Finder: %d cols removed, %d remaining\n", df_removed, ncols);
        }
    };

    {
        const int cols_before = ncols;
        // Step 1: Cheap reductions (cost + budget)
        apply_cost_reduction(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                             best_obj, tol, verbosity);
        apply_budget_pruning_preprocess(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                                        best_obj, tol, config.preprocess_time_limit, verbosity);

        // Step 2: Greedy multi-column dominance (fast, removes bulk of dominated cols)
        run_dominance_finder();

        // Step 3: Exact pairwise/triplet dominance on reduced set
        apply_dominance_reduction(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                                  "cost_driven", tol, config.preprocess_time_limit, verbosity);
        apply_dominance_reduction(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                                  config.preprocess_rules, tol, config.preprocess_time_limit, verbosity);

        // Step 4: Row reduction (essential columns, row domination, probing)
        auto rr = row_reduce(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                             tol, config.preprocess_time_limit, verbosity);
        if (rr.cols_fixed > 0 || rr.rows_removed > 0) {
            fixed_preprocess_cost += rr.fixed_cost;
            for (int c : rr.fixed_original_cols) fixed_original_cols.push_back(c);
            best_obj -= rr.fixed_cost;
        }
        // Iterate if row reduction changed the model
        if (rr.cols_fixed > 0 || rr.rows_removed > 0) {
            const int max_extra_rounds = 9;
            for (int round = 0; round < max_extra_rounds; ++round) {
                const int ncols_start = ncols;
                const int nrows_start = nrows;
                apply_cost_reduction(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                                     best_obj, tol, verbosity);
                apply_dominance_reduction(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                                          config.preprocess_rules, tol, config.preprocess_time_limit, verbosity);
                auto rr2 = row_reduce(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                                      tol, config.preprocess_time_limit, verbosity);
                if (rr2.cols_fixed > 0 || rr2.rows_removed > 0) {
                    fixed_preprocess_cost += rr2.fixed_cost;
                    for (int c : rr2.fixed_original_cols) fixed_original_cols.push_back(c);
                    best_obj -= rr2.fixed_cost;
                }
                if (ncols == ncols_start && nrows == nrows_start) break;
                if (verbosity >= 3)
                    fprintf(stderr, "  Preprocess round %d: %d/%d -> %d/%d (cols/rows)\n",
                            round + 2, ncols_start, nrows_start, ncols, nrows);
            }
        }
        if (ncols < cols_before && verbosity >= 2)
            fprintf(stderr, "  Pre-LP reduction: cols %d -> %d, rows %d\n", cols_before, ncols, nrows);
    }

    // Check if all rows are covered by essential columns
    if (nrows == 0) {
        const auto end_time = Clock::now();
        result.wall_time = std::chrono::duration<double>(end_time - start_time).count();
        result.primal_obj = fixed_preprocess_cost;
        result.dual_obj = fixed_preprocess_cost;
        result.mip_gap = 0.0;
        result.nodes_processed = 0;
        result.lp_solves = 0;
        result.status = "Optimal";
        result.solution.assign(static_cast<size_t>(ncols_input), 0.0);
        for (int c : fixed_original_cols) {
            if (c >= 0 && c < ncols_input)
                result.solution[static_cast<size_t>(c)] = 1.0;
        }
        if (verbosity >= 1)
            fprintf(stderr, "All rows covered by essential columns (cost %.12g)\n",
                    fixed_preprocess_cost);
        return result;
    }

    // ================================================================
    // Phase 3: Root LP solve
    // ================================================================
    if (verbosity >= 2) fprintf(stderr, "Phase 3: Root LP relaxation\n");

    BaseRelaxationModel base;
    auto build_base = [&]() {
        base.nrows = nrows;
        base.ncols = ncols;
        base.ncols_input = ncols_input;
        base.nnz = static_cast<int>(csr_vals.size());
        base.csr_inds = csr_inds;
        base.csr_offs = csr_offs;
        base.csr_vals = csr_vals;
        base.obj = obj;
        base.rhs.assign(static_cast<size_t>(nrows), 1.0);
        base.active_to_original = active_to_input;
        base.base_cuts.clear();
        base.build_transpose();
    };
    build_base();

    LpSolver lp;
    lp.build_model(base);
    int total_lp_solves = 0;

    // Root LP: use IPM for large models, simplex otherwise
    bool use_ipm_root = (config.root_lp_solver == "ipm") ||
        (config.root_lp_solver == "auto" && base.ncols >= 500 && base.nrows >= 200);
    LpSolution root_sol;
    if (use_ipm_root) {
        if (verbosity >= 2)
            fprintf(stderr, "  Root LP solver: IPM (barrier) with crossover\n");
        root_sol = lp.solve_ipm();
    } else {
        root_sol = lp.solve();
    }
    ++total_lp_solves;

    if (root_sol.solved && root_sol.optimal) {
        lp.save_basis();

        // Try heuristics on root LP
        BranchNodeState root_node;
        for (const auto &h : heuristics) {
            auto hr = h->tryBuild(root_sol.col_value, root_sol.row_dual, base, root_node, integ_tol);
            if (hr.feasible && hr.objective < best_obj - tol) {
                best_obj = hr.objective;
                adopt_incumbent_solution(best_solution, hr.solution, base.ncols, ncols_input, base.active_to_original);
                incumbent_source = std::string("root_") + hr.name;
            }
        }

        // Check if root LP is integral
        if (static_cast<int>(root_sol.col_value.size()) >= base.ncols &&
            is_binary_integral(root_sol.col_value, base.ncols, integ_tol) &&
            root_sol.primal_obj < best_obj - tol) {
            best_obj = root_sol.primal_obj;
            adopt_incumbent_solution(best_solution, root_sol.col_value, base.ncols, ncols_input, base.active_to_original);
            incumbent_source = "root_lp_exact";
        }

        double root_dual = root_sol.dual_obj;
        global_dual_bound_raw = root_dual;
        if (obj_is_integral)
            root_dual = tighten_dual_bound(root_dual, integ_tol);
        global_dual_bound = root_dual;

        if (verbosity >= 2)
            fprintf(stderr, "  Root LP: primal=%.12g dual=%.12g\n",
                    root_sol.primal_obj, root_sol.dual_obj);
    } else {
        if (verbosity >= 1)
            fprintf(stderr, "  Root LP did not converge\n");
    }

    // ================================================================
    // Phase 3.5: Lagrangian relaxation (subgradient optimization)
    // ================================================================
    if (root_sol.solved && root_sol.optimal) {
        const int lagr_max_iters = 500;
        const double lagr_time_limit = 5.0;
        auto lagr = lagrangian_relaxation(base, best_obj, root_sol.row_dual,
                                          lagr_max_iters, lagr_time_limit, verbosity);

        if (lagr.best_heuristic_obj < best_obj - tol &&
            !lagr.best_solution.empty()) {
            best_obj = lagr.best_heuristic_obj;
            adopt_incumbent_solution(best_solution, lagr.best_solution,
                                     base.ncols, ncols_input, base.active_to_original);
            incumbent_source = "lagrangian";
        }
    }

    // ================================================================
    // Phase 4-5: Post-LP reduction
    // ================================================================
    {
        const int cols_before = ncols;
        // Update working arrays from base
        csr_inds = base.csr_inds; csr_offs = base.csr_offs; csr_vals = base.csr_vals;
        obj = base.obj; active_to_input = base.active_to_original;
        ncols = base.ncols;

        // Iterative reduced-cost fixing: fix-resolve-fix cycles
        if (root_sol.solved && root_sol.optimal && std::isfinite(best_obj)) {
            const int max_rc_iters = 10;
            int total_rc_removed = 0;
            for (int rc_iter = 0; rc_iter < max_rc_iters; ++rc_iter) {
                const double gap = best_obj - root_sol.dual_obj;
                if (gap <= tol) break;

                auto rcosts = compute_reduced_costs(obj, root_sol.row_dual,
                    base, std::min(ncols, static_cast<int>(root_sol.col_value.size())));
                int rc_removed = 0;
                std::vector<char> rc_active(static_cast<size_t>(ncols), 1);
                for (int j = 0; j < ncols && j < static_cast<int>(rcosts.size()); ++j) {
                    if (root_sol.col_value[static_cast<size_t>(j)] < tol &&
                        rcosts[static_cast<size_t>(j)] > gap + tol) {
                        rc_active[static_cast<size_t>(j)] = 0;
                        ++rc_removed;
                    }
                }
                if (rc_removed == 0) break;

                // Rebuild without fixed columns
                std::vector<int> old_to_new(static_cast<size_t>(ncols), -1);
                std::vector<int> new_active_input;
                std::vector<double> new_obj_rc;
                int nc = 0;
                for (int j = 0; j < ncols; ++j) {
                    if (rc_active[static_cast<size_t>(j)]) {
                        old_to_new[static_cast<size_t>(j)] = nc;
                        new_active_input.push_back(active_to_input[static_cast<size_t>(j)]);
                        new_obj_rc.push_back(obj[static_cast<size_t>(j)]);
                        ++nc;
                    }
                }
                std::vector<int> new_inds_rc;
                std::vector<int> new_offs_rc;
                std::vector<double> new_vals_rc;
                new_offs_rc.push_back(0);
                for (int i = 0; i < nrows; ++i) {
                    for (int k = csr_offs[static_cast<size_t>(i)];
                         k < csr_offs[static_cast<size_t>(i) + 1]; ++k) {
                        const int c = csr_inds[static_cast<size_t>(k)];
                        if (c >= 0 && c < ncols) {
                            const int m = old_to_new[static_cast<size_t>(c)];
                            if (m >= 0) {
                                new_inds_rc.push_back(m);
                                new_vals_rc.push_back(csr_vals[static_cast<size_t>(k)]);
                            }
                        }
                    }
                    new_offs_rc.push_back(static_cast<int>(new_vals_rc.size()));
                }
                if (verbosity >= 3)
                    fprintf(stderr, "  RC fixing iter %d: %d cols removed, %d remaining\n",
                            rc_iter + 1, rc_removed, nc);
                ncols = nc;
                csr_inds = std::move(new_inds_rc);
                csr_offs = std::move(new_offs_rc);
                csr_vals = std::move(new_vals_rc);
                obj = std::move(new_obj_rc);
                active_to_input = std::move(new_active_input);
                total_rc_removed += rc_removed;

                // Rebuild base model and re-solve LP
                build_base();
                lp.rebuild_model(base);
                root_sol = lp.solve();
                ++total_lp_solves;
                if (!root_sol.solved || !root_sol.optimal) break;
                lp.save_basis();

                // Run heuristics on tighter LP (may improve incumbent → smaller gap)
                BranchNodeState rc_node;
                for (const auto &h : heuristics) {
                    auto hr = h->tryBuild(root_sol.col_value, root_sol.row_dual, base, rc_node, integ_tol);
                    if (hr.feasible && hr.objective < best_obj - tol) {
                        best_obj = hr.objective;
                        adopt_incumbent_solution(best_solution, hr.solution, base.ncols, ncols_input, base.active_to_original);
                        incumbent_source = std::string("rcfix_") + hr.name;
                    }
                }

                // Update working arrays from base for next iteration
                csr_inds = base.csr_inds; csr_offs = base.csr_offs; csr_vals = base.csr_vals;
                obj = base.obj; active_to_input = base.active_to_original;
                ncols = base.ncols;
            }
            if (total_rc_removed > 0 && verbosity >= 2)
                fprintf(stderr, "  Iterative RC fixing: %d cols removed total\n", total_rc_removed);
        }

        // Re-validate incumbent against current model after RC fixing.
        // RC fixing may have removed columns that were in the incumbent,
        // making the incumbent infeasible on the reduced model.
        if (std::isfinite(best_obj) && !best_solution.empty()) {
            // Check feasibility: every row must be covered by active incumbent columns
            bool incumbent_valid = true;
            double validated_obj = 0.0;
            for (int j = 0; j < ncols; ++j) {
                const int orig = active_to_input[static_cast<size_t>(j)];
                if (orig >= 0 && orig < ncols_input &&
                    best_solution[static_cast<size_t>(orig)] > 0.5)
                    validated_obj += obj[static_cast<size_t>(j)];
            }
            for (int i = 0; i < nrows && incumbent_valid; ++i) {
                double row_sum = 0.0;
                for (int k = csr_offs[static_cast<size_t>(i)];
                     k < csr_offs[static_cast<size_t>(i) + 1]; ++k) {
                    const int c = csr_inds[static_cast<size_t>(k)];
                    if (c >= 0 && c < ncols) {
                        const int orig = active_to_input[static_cast<size_t>(c)];
                        if (orig >= 0 && orig < ncols_input &&
                            best_solution[static_cast<size_t>(orig)] > 0.5)
                            row_sum += csr_vals[static_cast<size_t>(k)];
                    }
                }
                if (row_sum < 1.0 - tol) incumbent_valid = false;
            }
            if (!incumbent_valid || validated_obj > best_obj + tol) {
                if (verbosity >= 3)
                    fprintf(stderr, "  Incumbent invalidated after RC fixing (valid=%d, obj %.6g -> %.6g)\n",
                            incumbent_valid ? 1 : 0, best_obj, validated_obj);
                if (incumbent_valid) {
                    best_obj = validated_obj;
                } else {
                    // Incumbent infeasible on reduced model; re-run heuristics
                    best_obj = std::numeric_limits<double>::infinity();
                    if (root_sol.solved && root_sol.optimal) {
                        BranchNodeState reval_node;
                        for (const auto &h : heuristics) {
                            auto hr = h->tryBuild(root_sol.col_value, root_sol.row_dual, base, reval_node, integ_tol);
                            if (hr.feasible && hr.objective < best_obj - tol) {
                                best_obj = hr.objective;
                                adopt_incumbent_solution(best_solution, hr.solution, base.ncols, ncols_input, base.active_to_original);
                                incumbent_source = std::string("reval_") + hr.name;
                            }
                        }
                    }
                    if (verbosity >= 3)
                        fprintf(stderr, "  Re-computed incumbent: %.6g (%s)\n", best_obj, incumbent_source.c_str());
                }
            }
        }

        // LP-based probing: fix x_j=0, check if infeasible or bound > incumbent
        if (root_sol.solved && root_sol.optimal && std::isfinite(best_obj) && ncols <= 800) {
            const auto probe_start = Clock::now();
            const double probe_time_limit = 10.0;

            // base and LP are already current from iterative RC fixing (or initial build)
            {

                // Only probe columns with x_j > 0 in LP (fixing x_j=0 when already at 0 is a no-op)
                // Sort by decreasing LP value (most impactful to force to 0)
                std::vector<int> probe_order;
                for (int j = 0; j < ncols; ++j) {
                    if (root_sol.col_value[static_cast<size_t>(j)] > tol)
                        probe_order.push_back(j);
                }
                std::sort(probe_order.begin(), probe_order.end(),
                          [&](int a, int b) { return root_sol.col_value[static_cast<size_t>(a)] > root_sol.col_value[static_cast<size_t>(b)]; });

                std::vector<char> probe_fix_one(static_cast<size_t>(ncols), 0);
                int probe_fixed = 0;
                HighsBasis probe_basis = lp.get_basis();

                for (int j : probe_order) {
                    const double elapsed = std::chrono::duration<double>(Clock::now() - probe_start).count();
                    if (elapsed > probe_time_limit) break;

                    lp.restore_base_state();
                    lp.set_basis(probe_basis);
                    BranchDecision fix_zero{j, 0};
                    lp.apply_decisions({fix_zero});
                    auto probe_sol = lp.solve();
                    ++total_lp_solves;

                    if (!probe_sol.solved || probe_sol.infeasible ||
                        (probe_sol.optimal && probe_sol.primal_obj > best_obj + tol)) {
                        probe_fix_one[static_cast<size_t>(j)] = 1;
                        ++probe_fixed;
                    }
                }

                // Restore LP state
                lp.restore_base_state();
                lp.set_basis(probe_basis);

                if (probe_fixed > 0) {
                    if (verbosity >= 2)
                        fprintf(stderr, "  LP probing: %d cols fixed to 1\n", probe_fixed);

                    double probe_fixed_cost = 0.0;
                    double probe_incumbent_cost = 0.0;
                    std::vector<char> row_covered(static_cast<size_t>(nrows), 0);
                    for (int j = 0; j < ncols; ++j) {
                        if (!probe_fix_one[static_cast<size_t>(j)]) continue;
                        probe_fixed_cost += obj[static_cast<size_t>(j)];
                        const int orig = active_to_input[static_cast<size_t>(j)];
                        fixed_original_cols.push_back(orig);
                        if (!best_solution.empty() && orig >= 0 && orig < ncols_input &&
                            best_solution[static_cast<size_t>(orig)] > 0.5)
                            probe_incumbent_cost += obj[static_cast<size_t>(j)];
                        for (const auto &entry : base.cols_to_rows[static_cast<size_t>(j)]) {
                            if (entry.row >= 0 && entry.row < nrows && entry.val > 0.0)
                                row_covered[static_cast<size_t>(entry.row)] = 1;
                        }
                    }

                    // Rebuild without probed columns and their covered rows
                    std::vector<int> old_col_to_new(static_cast<size_t>(ncols), -1);
                    int new_ncols = 0;
                    std::vector<int> new_active;
                    std::vector<double> new_obj;
                    for (int j = 0; j < ncols; ++j) {
                        if (!probe_fix_one[static_cast<size_t>(j)]) {
                            old_col_to_new[static_cast<size_t>(j)] = new_ncols++;
                            new_active.push_back(active_to_input[static_cast<size_t>(j)]);
                            new_obj.push_back(obj[static_cast<size_t>(j)]);
                        }
                    }
                    int new_nrows = 0;
                    std::vector<int> new_csr_offs, new_csr_inds;
                    std::vector<double> new_csr_vals;
                    new_csr_offs.push_back(0);
                    for (int i = 0; i < nrows; ++i) {
                        if (row_covered[static_cast<size_t>(i)]) continue;
                        for (int k = csr_offs[static_cast<size_t>(i)]; k < csr_offs[static_cast<size_t>(i) + 1]; ++k) {
                            const int c = csr_inds[static_cast<size_t>(k)];
                            if (c >= 0 && c < ncols) {
                                const int m = old_col_to_new[static_cast<size_t>(c)];
                                if (m >= 0) {
                                    new_csr_inds.push_back(m);
                                    new_csr_vals.push_back(csr_vals[static_cast<size_t>(k)]);
                                }
                            }
                        }
                        new_csr_offs.push_back(static_cast<int>(new_csr_vals.size()));
                        ++new_nrows;
                    }

                    ncols = new_ncols;
                    nrows = new_nrows;
                    csr_inds = std::move(new_csr_inds);
                    csr_offs = std::move(new_csr_offs);
                    csr_vals = std::move(new_csr_vals);
                    obj = std::move(new_obj);
                    active_to_input = std::move(new_active);

                    fixed_preprocess_cost += probe_fixed_cost;
                    best_obj -= probe_incumbent_cost;
                    if (std::isfinite(global_dual_bound))
                        global_dual_bound -= probe_fixed_cost;
                    if (std::isfinite(global_dual_bound_raw))
                        global_dual_bound_raw -= probe_fixed_cost;

                    // Rebuild and re-solve for heuristics on reduced model
                    build_base();
                    lp.rebuild_model(base);
                    root_sol = lp.solve();
                    ++total_lp_solves;
                    if (root_sol.solved && root_sol.optimal) {
                        lp.save_basis();
                        BranchNodeState probe_node;
                        for (const auto &h : heuristics) {
                            auto hr = h->tryBuild(root_sol.col_value, root_sol.row_dual, base, probe_node, integ_tol);
                            if (hr.feasible && hr.objective < best_obj - tol) {
                                best_obj = hr.objective;
                                adopt_incumbent_solution(best_solution, hr.solution, base.ncols, ncols_input, base.active_to_original);
                                incumbent_source = std::string("probe_") + hr.name;
                            }
                        }
                        csr_inds = base.csr_inds; csr_offs = base.csr_offs; csr_vals = base.csr_vals;
                        obj = base.obj; active_to_input = base.active_to_original;
                        ncols = base.ncols;
                    }
                }
            }
        }

        // Post-LP reduction (first round)
        apply_cost_reduction(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                             best_obj, tol, verbosity);
        apply_budget_pruning_preprocess(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                                        best_obj, tol, config.preprocess_time_limit, verbosity);
        apply_dominance_reduction(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                                  config.preprocess_rules, tol, config.preprocess_time_limit, verbosity);
        {
            auto rr2 = row_reduce(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                                  tol, config.preprocess_time_limit, verbosity);
            if (rr2.cols_fixed > 0 || rr2.rows_removed > 0) {
                fixed_preprocess_cost += rr2.fixed_cost;
                for (int c : rr2.fixed_original_cols) fixed_original_cols.push_back(c);
                best_obj -= rr2.fixed_cost;
                if (std::isfinite(global_dual_bound))
                    global_dual_bound -= rr2.fixed_cost;
                if (std::isfinite(global_dual_bound_raw))
                    global_dual_bound_raw -= rr2.fixed_cost;
            }
            // Iterate only if row reduction found something
            if (rr2.cols_fixed > 0 || rr2.rows_removed > 0) {
                const int max_extra_rounds = 9;
                for (int round = 0; round < max_extra_rounds; ++round) {
                    const int ncols_start = ncols;
                    const int nrows_start = nrows;
                    apply_cost_reduction(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                                         best_obj, tol, verbosity);
                    apply_dominance_reduction(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                                              config.preprocess_rules, tol, config.preprocess_time_limit, verbosity);
                    auto rr3 = row_reduce(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                                          tol, config.preprocess_time_limit, verbosity);
                    if (rr3.cols_fixed > 0 || rr3.rows_removed > 0) {
                        fixed_preprocess_cost += rr3.fixed_cost;
                        for (int c : rr3.fixed_original_cols) fixed_original_cols.push_back(c);
                        best_obj -= rr3.fixed_cost;
                        if (std::isfinite(global_dual_bound))
                            global_dual_bound -= rr3.fixed_cost;
                        if (std::isfinite(global_dual_bound_raw))
                            global_dual_bound_raw -= rr3.fixed_cost;
                    }
                    if (ncols == ncols_start && nrows == nrows_start) break;
                    if (verbosity >= 3)
                        fprintf(stderr, "  Post-LP round %d: %d/%d -> %d/%d (cols/rows)\n",
                                round + 2, ncols_start, nrows_start, ncols, nrows);
                }
            }
        }
        if (ncols < cols_before && verbosity >= 2)
            fprintf(stderr, "  Post-LP reduction: cols %d -> %d, rows %d\n", cols_before, ncols, nrows);
    }

    // ================================================================
    // Phase 6: Build base model for BnB
    // ================================================================
    build_base();
    lp.rebuild_model(base);

    if (root_sol.solved && root_sol.optimal) {
        // Re-solve to get warm basis for reduced model
        root_sol = lp.solve();
        ++total_lp_solves;
        if (root_sol.solved && root_sol.optimal) lp.save_basis();
    }

    if (verbosity >= 2) {
        const double density = (base.nrows > 0 && base.ncols > 0)
            ? 100.0 * base.nnz / (static_cast<double>(base.nrows) * base.ncols)
            : 0.0;
        fprintf(stderr, "BnB base model: %d rows x %d cols, nnz=%d (%.1f%% dense)\n",
                base.nrows, base.ncols, base.nnz, density);
    }

    // ================================================================
    // Phase 6.5: Root cut separation rounds
    // ================================================================
    int root_cuts_added = 0;
    if (config.cuts_enabled && config.cut_rounds_root > 0 && root_sol.solved && root_sol.optimal) {
        auto separators = make_cut_separators();
        LpSolution cut_sol = root_sol;

        // Save pre-cut model state to restore if cuts don't help
        const double pre_cut_tightened_dual = global_dual_bound;
        const int pre_cut_nrows = base.nrows;
        const int pre_cut_nnz = base.nnz;
        const auto pre_cut_csr_inds = base.csr_inds;
        const auto pre_cut_csr_offs = base.csr_offs;
        const auto pre_cut_csr_vals = base.csr_vals;
        const auto pre_cut_rhs = base.rhs;
        const auto pre_cut_base_cuts = base.base_cuts;

        for (int round = 0; round < config.cut_rounds_root; ++round) {
            // Check integrality
            if (static_cast<int>(cut_sol.col_value.size()) >= base.ncols &&
                is_binary_integral(cut_sol.col_value, base.ncols, integ_tol) &&
                cut_sol.primal_obj < best_obj - tol) {
                best_obj = cut_sol.primal_obj;
                adopt_incumbent_solution(best_solution, cut_sol.col_value, base.ncols, ncols_input, base.active_to_original);
                incumbent_source = "cut_round_exact";
                if (verbosity >= 2)
                    fprintf(stderr, "  Cut round %d: LP integral, incumbent %.12g\n", round + 1, best_obj);
                break;
            }

            // Heuristics
            BranchNodeState empty_branch;
            for (const auto &h : heuristics) {
                auto hr = h->tryBuild(cut_sol.col_value, cut_sol.row_dual, base, empty_branch, integ_tol);
                if (hr.feasible && hr.objective < best_obj - tol) {
                    best_obj = hr.objective;
                    adopt_incumbent_solution(best_solution, hr.solution, base.ncols, ncols_input, base.active_to_original);
                    incumbent_source = std::string("cut_") + hr.name;
                }
            }

            // Update dual bound
            if (cut_sol.optimal && std::isfinite(cut_sol.dual_obj)) {
                double cut_dual = cut_sol.dual_obj;
                if (cut_dual > global_dual_bound_raw) global_dual_bound_raw = cut_dual;
                if (obj_is_integral) cut_dual = tighten_dual_bound(cut_dual, integ_tol);
                if (cut_dual > global_dual_bound) global_dual_bound = cut_dual;
            }

            // Separate
            std::vector<CutConstraint> round_cuts;
            for (const auto &sep : separators) {
                auto sep_cuts = sep->separate(cut_sol.col_value, cut_sol.row_dual,
                                              base, base.ncols, integ_tol, best_obj);
                for (auto &c : sep_cuts) {
                    if (static_cast<int>(round_cuts.size()) >= config.max_cuts_per_round) break;
                    round_cuts.push_back(std::move(c));
                }
                if (static_cast<int>(round_cuts.size()) >= config.max_cuts_per_round) break;
            }
            if (round_cuts.empty()) {
                if (verbosity >= 3)
                    fprintf(stderr, "  Cut round %d: no violated cuts\n", round + 1);
                break;
            }

            append_cuts_to_base_model(base, round_cuts);
            root_cuts_added += static_cast<int>(round_cuts.size());
            if (verbosity >= 3)
                fprintf(stderr, "  Cut round %d: added %d cuts (model %d rows)\n",
                        round + 1, static_cast<int>(round_cuts.size()), base.nrows);

            // Rebuild LP with cuts
            lp.rebuild_model(base);
            cut_sol = lp.solve();
            ++total_lp_solves;
            if (!cut_sol.solved || !cut_sol.optimal) break;
            lp.save_basis();
        }

        // Check if cuts actually improved the tightened dual bound.
        // For integer objectives, cuts that don't push the LP past the next integer
        // provide no benefit and can disrupt branching topology.
        if (root_cuts_added > 0 && global_dual_bound <= pre_cut_tightened_dual + tol) {
            // Cuts didn't improve the tightened bound — restore pre-cut model
            base.nrows = pre_cut_nrows;
            base.nnz = pre_cut_nnz;
            base.csr_inds = pre_cut_csr_inds;
            base.csr_offs = pre_cut_csr_offs;
            base.csr_vals = pre_cut_csr_vals;
            base.rhs = pre_cut_rhs;
            base.base_cuts = pre_cut_base_cuts;
            lp.rebuild_model(base);
            root_cuts_added = 0;
            if (verbosity >= 2)
                fprintf(stderr, "  Root cuts undone: no tightened dual improvement\n");
        } else if (root_cuts_added > 0 && verbosity >= 2) {
            fprintf(stderr, "  Root cuts: %d total (dual %.12g -> %.12g)\n",
                    root_cuts_added, pre_cut_tightened_dual, global_dual_bound);
        }
    }

    // ================================================================
    // Phase 6.6: Root Balas cover cuts
    // ================================================================
    int root_balas_cuts_added = 0;
    if (config.balas_enabled && root_sol.solved && root_sol.optimal) {
        // Reuse existing LP solution (avoid unnecessary re-solve that perturbs basis)
        LpSolution balas_sol = lp.solve();
        ++total_lp_solves;
        if (balas_sol.solved && balas_sol.optimal) {
            auto rcosts = compute_reduced_costs(base.obj, balas_sol.row_dual, base, base.ncols);
            auto br = balas_branch_generate(balas_sol.col_value, balas_sol.row_dual, rcosts,
                                            best_obj, base, base.ncols, config.balas_max_branches,
                                            integ_tol, false);
            if (!br.sets.empty()) {
                std::vector<CutConstraint> balas_cuts;
                for (const auto &R_k : br.sets) {
                    // Only add cover cuts that are violated by the LP
                    double lhs_val = 0.0;
                    for (const int j : R_k) {
                        if (j >= 0 && j < base.ncols)
                            lhs_val += balas_sol.col_value[static_cast<size_t>(j)];
                    }
                    if (lhs_val >= 1.0 - integ_tol) continue;

                    CutConstraint cut;
                    cut.indices.reserve(R_k.size());
                    cut.values.reserve(R_k.size());
                    for (const int j : R_k) {
                        cut.indices.push_back(j);
                        cut.values.push_back(1.0);
                    }
                    cut.rhs = 1.0;
                    cut.type = ">=";
                    balas_cuts.push_back(std::move(cut));
                }
                if (!balas_cuts.empty()) {
                    append_cuts_to_base_model(base, balas_cuts);
                    root_balas_cuts_added = static_cast<int>(balas_cuts.size());
                    lp.rebuild_model(base);
                    if (verbosity >= 2)
                        fprintf(stderr, "  Root Balas cover cuts: %d added (from %d sets)\n",
                                root_balas_cuts_added, static_cast<int>(br.sets.size()));
                }
            }
        }
    }

    // ================================================================
    // Phase 6.7: Post-cut budget pruning
    // ================================================================
    if (std::isfinite(best_obj)) {
        ModelReductionResult budget_red = reduce_base_model_budget_pruning(
            base, best_obj, tol, config.preprocess_time_limit);
        if (budget_red.columns_removed > 0) {
            if (verbosity >= 3)
                fprintf(stderr, "  Post-cut budget pruning: %d cols removed\n",
                        budget_red.columns_removed);
            lp.rebuild_model(base);
        }
    }

    // ================================================================
    // BnB main loop
    // ================================================================
    const bool use_reliability = (config.branch_strategy == "reliability");
    auto selector = make_branch_selector(
        use_reliability ? "most_fractional" : config.branch_strategy);
    PseudocostState pc_state;
    pc_state.init(base.ncols);

    // Cache adjacency for reliability branching propagation
    ScpAdjacency cached_adj = build_adjacency(base);
    bool adj_valid = true;

    // Build orbital fixing state for symmetry breaking
    OrbitalFixingState orbital_state = build_orbital_state(base);
    bool orbital_valid = true;
    if (verbosity >= 2 && !orbital_state.cost_groups.empty()) {
        int total_cols_in_groups = 0;
        for (const auto &g : orbital_state.cost_groups)
            total_cols_in_groups += static_cast<int>(g.size());
        fprintf(stderr, "  Symmetry: %d cost groups, %d cols in groups\n",
                static_cast<int>(orbital_state.cost_groups.size()), total_cols_in_groups);
    }

    // Decomposition state (computed after all root preprocessing so linking
    // columns are only among truly free variables, not RC-fixed ones)
    DecompositionState decomp_state;
    bool decomp_dynamic_eligible = (config.decomposition_mode == "auto");
    if (config.decomposition_mode == "force" ||
        (config.decomposition_mode == "auto" && base.ncols > 1000)) {
        decomp_state = compute_decomposition(base, config.decomposition_max_linkers, verbosity);
        decomp_dynamic_eligible = false;
    }

    std::vector<BranchNodeState> nodes;
    std::unordered_map<int, HighsBasis> node_bases;
    {
        BranchNodeState root_state;
        root_state.parent_dual_bound = std::isfinite(global_dual_bound)
                                           ? global_dual_bound
                                           : -std::numeric_limits<double>::infinity();
        root_state.parent_dual_bound_raw = std::isfinite(global_dual_bound_raw)
                                               ? global_dual_bound_raw
                                               : -std::numeric_limits<double>::infinity();
        nodes.push_back(root_state);
    }

    BestBoundFrontier frontier;
    frontier.init(&nodes);
    frontier.push(0);

    int processed_nodes = 0;
    bool gap_tolerance_reached = false;
    bool hard_time_limit_reached = false;
    bool frontier_exhausted = false;

    const int gap_stagnation_window = config.gap_stagnation_window;
    double best_mip_gap_seen = std::numeric_limits<double>::infinity();
    int node_at_last_gap_improvement = 0;

    double cut_accumulator = 0.0;
    double balas_accumulator = 0.0;
    double lagrangian_accumulator = 0.0;
    double diving_accumulator = 0.0;
    double rins_accumulator = 0.0;
    // Incumbent in active column space (for RINS comparison with LP)
    std::vector<double> incumbent_active;
    size_t frontier_at_last_stagnation = 0;
    int frontier_shrink_streak = 0;
    bool proving_phase = false; // true when frontier monotonically decreasing
    bool force_aggressive_branching = false;
    auto cut_separators = make_cut_separators();

    const auto bnb_start = Clock::now();
    double next_log_sec = config.log_interval_seconds;
    int log_event_count = 0;

    if (verbosity >= 2)
        fprintf(stderr, "Branch-and-bound started (max_nodes=%d)\n", config.max_nodes);

    while (processed_nodes < config.max_nodes) {
        // Time limit check
        const auto now = Clock::now();
        const double elapsed = std::chrono::duration<double>(now - start_time).count();
        if (config.time_limit_seconds > 0.0 && elapsed >= config.time_limit_seconds) {
            hard_time_limit_reached = true;
            if (verbosity >= 1)
                fprintf(stderr, "  [%8.3fs] Time limit reached\n", elapsed);
            break;
        }

        // Gap check
        if (std::isfinite(best_obj) && std::isfinite(global_dual_bound)) {
            const double gap = compute_mip_gap(best_obj, global_dual_bound);
            if (std::isfinite(gap) && gap <= config.mip_gap_tol) {
                gap_tolerance_reached = true;
                if (verbosity >= 1) {
                    const double t = std::chrono::duration<double>(now - start_time).count();
                    fprintf(stderr, "  [%8.3fs] MIP gap %.8f%% within tolerance; optimal\n", t, gap * 100.0);
                }
                break;
            }
        }

        // Log
        if (config.log_interval_seconds > 0.0) {
            const double bnb_elapsed = std::chrono::duration<double>(now - bnb_start).count();
            if (bnb_elapsed >= next_log_sec) {
                if (verbosity >= 2) {
                    const double gap = compute_mip_gap(best_obj, global_dual_bound_raw);
                    fprintf(stderr, "  [%8.3fs] nodes=%6d frontier=%6zu lp=%6d incumbent=%.8f dual=%.8f gap=",
                            elapsed, processed_nodes, frontier.size(), total_lp_solves,
                            best_obj, global_dual_bound_raw);
                    if (std::isfinite(gap)) fprintf(stderr, "%.8f%%\n", gap * 100.0);
                    else fprintf(stderr, "inf\n");
                    ++log_event_count;
                    if (log_event_count % 10 == 0) {
                        int total_cuts = base.nrows - nrows;
                        fprintf(stderr, "           model: %d rows x %d cols, %d cuts\n",
                                base.nrows, base.ncols, total_cuts);
                    }
                }
                next_log_sec = std::chrono::duration<double>(now - bnb_start).count() + config.log_interval_seconds;
            }
        }

        // Get next node
        if (frontier.empty()) {
            frontier_exhausted = true;
            break;
        }
        const int node_id = frontier.pop();
        const BranchNodeState branch_node = nodes[static_cast<size_t>(node_id)];

        // Reclaim memory from processed node
        {
            auto &nd = nodes[static_cast<size_t>(node_id)];
            nd.decisions.clear(); nd.decisions.shrink_to_fit();
            nd.cuts.clear(); nd.cuts.shrink_to_fit();
        }

        // Retrieve per-node basis if available
        auto basis_it = node_bases.find(node_id);
        const bool has_node_basis = (basis_it != node_bases.end());

        // Bound check
        if (branch_node.parent_dual_bound >= best_obj - tol) {
            if (has_node_basis) node_bases.erase(basis_it);
            continue;
        }

        // Node-level propagation: detect infeasibility and find implied fixings.
        // Only worthwhile when enough variables are fixed to create essential columns.
        const bool do_propagation = (static_cast<int>(branch_node.decisions.size()) >= 3);
        std::vector<BranchDecision> augmented_decisions;
        if (do_propagation) {
            if (!adj_valid) {
                cached_adj = build_adjacency(base);
                adj_valid = true;
            }
            auto prop = propagate_fixings(cached_adj, branch_node.decisions, -1, 0);
            if (prop.infeasible) {
                ++processed_nodes;
                continue;
            }
            if (!prop.implied.empty()) {
                augmented_decisions = branch_node.decisions;
                for (const auto &imp : prop.implied)
                    augmented_decisions.push_back(imp);
            }
        }

        // Orbital fixing: fix equivalent columns in the same symmetry class
        if (!orbital_state.disabled &&
            static_cast<int>(branch_node.decisions.size()) >= 2) {
            if (!orbital_valid) {
                orbital_state = build_orbital_state(base);
                orbital_valid = true;
            }
            if (!orbital_state.cost_groups.empty()) {
                if (!adj_valid) {
                    cached_adj = build_adjacency(base);
                    adj_valid = true;
                }
                const auto &effective = augmented_decisions.empty()
                                            ? branch_node.decisions
                                            : augmented_decisions;
                auto orbit_fix = compute_orbital_fixings(cached_adj, effective, orbital_state);
                if (!orbit_fix.empty()) {
                    if (augmented_decisions.empty())
                        augmented_decisions = branch_node.decisions;
                    for (const auto &f : orbit_fix)
                        augmented_decisions.push_back(f);
                    orbital_state.total_fixed += static_cast<int>(orbit_fix.size());
                }
                // Auto-disable if no fixings found after initial probe period
                ++orbital_state.calls;
                if (orbital_state.calls >= 100 && orbital_state.total_fixed == 0)
                    orbital_state.disabled = true;
            }
        }

        // Decomposition: if all linking columns are decided, solve sub-problems.
        // Pass augmented_decisions (branch + propagated + orbital + RC) to
        // solve_decomposed so sub-problems see the full node state.
        if (decomp_state.enabled) {
            const auto &eff_dec = augmented_decisions.empty()
                                      ? branch_node.decisions
                                      : augmented_decisions;
            if (all_linkers_fixed(decomp_state, eff_dec)) {
                const double bnb_elapsed = std::chrono::duration<double>(
                    Clock::now() - start_time).count();
                const double remaining = (config.time_limit_seconds > 0)
                    ? std::max(1.0, config.time_limit_seconds - bnb_elapsed)
                    : 60.0;

                auto dr = solve_decomposed(base, eff_dec, decomp_state,
                                           best_obj, config, remaining, verbosity);
                if (dr.solved) {
                    // Sub-problems fully explored this subtree
                    if (dr.improved) {
                        best_obj = dr.combined_obj;
                        adopt_incumbent_solution(best_solution, dr.combined_solution,
                                                 base.ncols, ncols_input,
                                                 base.active_to_original);
                        incumbent_source = "decomposition";
                        if (verbosity >= 2) {
                            const double t = std::chrono::duration<double>(
                                Clock::now() - start_time).count();
                            fprintf(stderr,
                                "* [%8.3fs] New incumbent: %.8f (from decomposition)\n",
                                t, best_obj);
                        }
                        node_at_last_gap_improvement = processed_nodes;
                        frontier.prune(best_obj, tol, verbosity, &node_bases);
                        if (mid_bnb_column_removal(base, best_obj, tol, frontier,
                                                    nodes, verbosity, &pc_state) > 0) {
                            lp.rebuild_model_keep_basis(base);
                            adj_valid = false; orbital_valid = false;
                        }
                        if (mid_bnb_budget_pruning(base, best_obj, tol,
                                                    config.preprocess_time_limit,
                                                    frontier, nodes, verbosity,
                                                    &pc_state) > 0) {
                            lp.rebuild_model_keep_basis(base);
                            adj_valid = false; orbital_valid = false;
                        }
                    }
                    continue; // skip LP solve/branch — subtree handled
                }
                // Decomposition failed (infeasible block or no split):
                // fall through to normal LP solve + branching
            }
        }

        // Solve LP with per-node basis warm-start
        lp.apply_decisions(augmented_decisions.empty() ? branch_node.decisions : augmented_decisions);
        lp.add_cuts(branch_node.cuts);
        if (has_node_basis) {
            lp.set_basis(basis_it->second);
            node_bases.erase(basis_it);
        } else {
            lp.restore_basis();
        }
        LpSolution sol = lp.solve();
        ++total_lp_solves;
        lp.save_basis();
        HighsBasis node_basis = lp.get_basis();
        lp.restore_base_state();

        if (!sol.solved || sol.infeasible) {
            ++processed_nodes;
            continue;
        }

        ++processed_nodes;

        double node_dual_bound_raw = sol.optimal ? sol.dual_obj : branch_node.parent_dual_bound_raw;
        double node_dual_bound = sol.optimal ? sol.dual_obj : branch_node.parent_dual_bound;
        if (obj_is_integral && sol.optimal && std::isfinite(node_dual_bound))
            node_dual_bound = tighten_dual_bound(node_dual_bound, integ_tol);

        const bool dual_improved = sol.optimal &&
            (node_dual_bound > branch_node.parent_dual_bound + tol);

        // Heuristics
        if (processed_nodes == 1 ||
            (config.heuristic_frequency > 0 && processed_nodes % config.heuristic_frequency == 0) ||
            dual_improved) {
            bool incumbent_improved = false;
            for (const auto &h : heuristics) {
                auto hr = h->tryBuild(sol.col_value, sol.row_dual, base, branch_node, integ_tol);
                if (hr.feasible && hr.objective < best_obj - tol) {
                    best_obj = hr.objective;
                    adopt_incumbent_solution(best_solution, hr.solution, base.ncols, ncols_input, base.active_to_original);
                    incumbent_source = hr.name;
                    incumbent_improved = true;
                    if (verbosity >= 2) {
                        const double t = std::chrono::duration<double>(Clock::now() - start_time).count();
                        fprintf(stderr, "* [%8.3fs] New incumbent: %.8f (from %s)\n", t, best_obj, hr.name.c_str());
                    }
                    break;
                }
            }
            if (incumbent_improved) {
                node_at_last_gap_improvement = processed_nodes;
                frontier.prune(best_obj, tol, verbosity, &node_bases);
                if (mid_bnb_column_removal(base, best_obj, tol, frontier, nodes, verbosity, &pc_state) > 0) {
                    lp.rebuild_model_keep_basis(base); adj_valid = false; orbital_valid = false;
                }
                if (mid_bnb_budget_pruning(base, best_obj, tol, config.preprocess_time_limit, frontier, nodes, verbosity, &pc_state) > 0) {
                    lp.rebuild_model_keep_basis(base); adj_valid = false; orbital_valid = false;
                }
            }
        }

        // Node-level Lagrangian heuristic
        if (config.lagrangian_frequency > 0.0 && sol.optimal &&
            static_cast<int>(sol.row_dual.size()) >= base.nrows) {
            lagrangian_accumulator += config.lagrangian_frequency;
            if (lagrangian_accumulator >= 1.0) {
                lagrangian_accumulator -= 1.0;

                const auto &eff_dec = augmented_decisions.empty()
                                          ? branch_node.decisions
                                          : augmented_decisions;

                auto lagr = lagrangian_relaxation_at_node(
                    base, best_obj, sol.row_dual, eff_dec,
                    50,   // max_iterations (reduced from 500)
                    0.5); // time_limit (reduced from 5.0)

                if (lagr.best_heuristic_obj < best_obj - tol &&
                    !lagr.best_solution.empty()) {
                    best_obj = lagr.best_heuristic_obj;
                    adopt_incumbent_solution(best_solution, lagr.best_solution,
                                             base.ncols, ncols_input,
                                             base.active_to_original);
                    incumbent_source = "node_lagrangian";
                    if (verbosity >= 2) {
                        const double t = std::chrono::duration<double>(
                            Clock::now() - start_time).count();
                        fprintf(stderr,
                            "* [%8.3fs] New incumbent: %.8f (from node_lagrangian)\n",
                            t, best_obj);
                    }
                    node_at_last_gap_improvement = processed_nodes;
                    frontier.prune(best_obj, tol, verbosity, &node_bases);
                    if (mid_bnb_column_removal(base, best_obj, tol, frontier,
                                                nodes, verbosity, &pc_state) > 0) {
                        lp.rebuild_model_keep_basis(base);
                        adj_valid = false; orbital_valid = false;
                    }
                    if (mid_bnb_budget_pruning(base, best_obj, tol,
                                                config.preprocess_time_limit,
                                                frontier, nodes, verbosity,
                                                &pc_state) > 0) {
                        lp.rebuild_model_keep_basis(base);
                        adj_valid = false; orbital_valid = false;
                    }
                }
            }
        }

        // Diving heuristic (skip in proving phase — frontier shrinking)
        if (config.diving_frequency > 0.0 && sol.optimal && !proving_phase &&
            node_dual_bound < best_obj - tol) {
            diving_accumulator += config.diving_frequency;
            if (diving_accumulator >= 1.0) {
                diving_accumulator -= 1.0;

                auto dive = run_diving_heuristics(
                    lp, base, branch_node, sol, best_obj, integ_tol,
                    config.diving_max_lp_solves);
                total_lp_solves += dive.lp_solves;

                if (dive.found && dive.objective < best_obj - tol) {
                    best_obj = dive.objective;
                    adopt_incumbent_solution(best_solution, dive.solution,
                                             base.ncols, ncols_input,
                                             base.active_to_original);
                    incumbent_source = dive.strategy_name;
                    if (verbosity >= 2) {
                        const double t = std::chrono::duration<double>(
                            Clock::now() - start_time).count();
                        fprintf(stderr,
                            "* [%8.3fs] New incumbent: %.8f (from %s)\n",
                            t, best_obj, dive.strategy_name.c_str());
                    }
                    node_at_last_gap_improvement = processed_nodes;
                    frontier.prune(best_obj, tol, verbosity, &node_bases);
                    if (mid_bnb_column_removal(base, best_obj, tol, frontier,
                                                nodes, verbosity, &pc_state) > 0) {
                        lp.rebuild_model_keep_basis(base);
                        adj_valid = false; orbital_valid = false;
                    }
                    if (mid_bnb_budget_pruning(base, best_obj, tol,
                                                config.preprocess_time_limit,
                                                frontier, nodes, verbosity,
                                                &pc_state) > 0) {
                        lp.rebuild_model_keep_basis(base);
                        adj_valid = false; orbital_valid = false;
                    }
                }
            }
        }

        // RINS heuristic (skip in proving phase — frontier shrinking)
        if (config.rins_frequency > 0.0 && sol.optimal && !proving_phase &&
            !best_solution.empty() && node_dual_bound < best_obj - tol) {
            rins_accumulator += config.rins_frequency;
            if (rins_accumulator >= 1.0) {
                rins_accumulator -= 1.0;

                // Build incumbent in active column space (lazy)
                incumbent_active.assign(static_cast<size_t>(base.ncols), 0.0);
                for (int j = 0; j < base.ncols; ++j) {
                    const int orig = base.active_to_original[static_cast<size_t>(j)];
                    if (orig >= 0 && orig < ncols_input &&
                        best_solution[static_cast<size_t>(orig)] > 0.5)
                        incumbent_active[static_cast<size_t>(j)] = 1.0;
                }

                const double bnb_elapsed = std::chrono::duration<double>(
                    Clock::now() - start_time).count();
                const double remaining = (config.time_limit_seconds > 0)
                    ? std::max(1.0, config.time_limit_seconds - bnb_elapsed)
                    : 30.0;

                auto rr = run_rins(base, sol, incumbent_active, best_obj,
                                   config, std::min(remaining, 30.0), verbosity);
                total_lp_solves += rr.sub_lp_solves;

                if (rr.found && rr.objective < best_obj - tol) {
                    best_obj = rr.objective;
                    adopt_incumbent_solution(best_solution, rr.solution,
                                             base.ncols, ncols_input,
                                             base.active_to_original);
                    incumbent_source = "rins";
                    if (verbosity >= 2) {
                        const double t = std::chrono::duration<double>(
                            Clock::now() - start_time).count();
                        fprintf(stderr,
                            "* [%8.3fs] New incumbent: %.8f (from rins)\n",
                            t, best_obj);
                    }
                    node_at_last_gap_improvement = processed_nodes;
                    frontier.prune(best_obj, tol, verbosity, &node_bases);
                    if (mid_bnb_column_removal(base, best_obj, tol, frontier,
                                                nodes, verbosity, &pc_state) > 0) {
                        lp.rebuild_model_keep_basis(base);
                        adj_valid = false; orbital_valid = false;
                    }
                    if (mid_bnb_budget_pruning(base, best_obj, tol,
                                                config.preprocess_time_limit,
                                                frontier, nodes, verbosity,
                                                &pc_state) > 0) {
                        lp.rebuild_model_keep_basis(base);
                        adj_valid = false; orbital_valid = false;
                    }
                }
            }
        }

        // Pruning by bound
        if (node_dual_bound >= best_obj - tol) continue;

        // Integrality check
        if (sol.optimal && is_binary_integral(sol.col_value, base.ncols, integ_tol)) {
            if (sol.primal_obj < best_obj - tol) {
                best_obj = sol.primal_obj;
                node_at_last_gap_improvement = processed_nodes;
                adopt_incumbent_solution(best_solution, sol.col_value, base.ncols, ncols_input, base.active_to_original);
                incumbent_source = "exact_node";
                if (verbosity >= 2) {
                    const double t = std::chrono::duration<double>(Clock::now() - start_time).count();
                    fprintf(stderr, "* [%8.3fs] New incumbent: %.8f (exact)\n", t, best_obj);
                }
                frontier.prune(best_obj, tol, verbosity, &node_bases);
                if (mid_bnb_column_removal(base, best_obj, tol, frontier, nodes, verbosity, &pc_state) > 0) {
                    lp.rebuild_model_keep_basis(base); adj_valid = false; orbital_valid = false;
                }
                if (mid_bnb_budget_pruning(base, best_obj, tol, config.preprocess_time_limit, frontier, nodes, verbosity, &pc_state) > 0) {
                    lp.rebuild_model_keep_basis(base); adj_valid = false; orbital_valid = false;
                }
            }
            continue;
        }

        // Compute reduced costs once per node (used by RC fixing and Balas)
        std::vector<double> rcosts;
        const bool need_rcosts = sol.optimal &&
            ((std::isfinite(best_obj) && node_dual_bound < best_obj - tol) ||
             (config.balas_enabled && force_aggressive_branching));
        if (need_rcosts)
            rcosts = compute_reduced_costs(base.obj, sol.row_dual, base, base.ncols);

        // Node-level reduced cost fixing: fix variables that cannot improve
        std::vector<BranchDecision> rc_fixings;
        if (!rcosts.empty() && std::isfinite(best_obj) && node_dual_bound < best_obj - tol) {
            const double node_gap = best_obj - node_dual_bound;
            for (int j = 0; j < base.ncols; ++j) {
                if (sol.col_value[static_cast<size_t>(j)] < integ_tol &&
                    rcosts[static_cast<size_t>(j)] > node_gap + tol) {
                    rc_fixings.push_back({j, 0});
                }
            }
        }

        // Branch
        auto fractional = collect_fractional_candidates(sol.col_value, base.ncols, integ_tol);
        if (fractional.empty()) continue;

        bool used_balas = false;
        if (config.balas_enabled && force_aggressive_branching && !rcosts.empty()) {
            auto br = balas_branch_generate(sol.col_value, sol.row_dual, rcosts,
                                            best_obj, base, base.ncols, config.balas_max_branches,
                                            integ_tol, true);
            if (br.use_br1) {
                auto children = balas_br1_children(branch_node, br.sets, node_dual_bound, node_dual_bound_raw);
                int enqueued = 0;
                for (auto &child : children) {
                    if (!is_node_provably_infeasible(child, base)) {
                        const int child_id = static_cast<int>(nodes.size());
                        nodes.push_back(std::move(child));
                        node_bases[child_id] = node_basis;
                        frontier.push(child_id);
                        ++enqueued;
                    }
                }
                used_balas = (enqueued > 0);
                if (verbosity >= 3) {
                    const double t = std::chrono::duration<double>(Clock::now() - start_time).count();
                    fprintf(stderr, "  [%8.3fs] Aggressive BR1: %d children from %d sets\n",
                            t, enqueued, static_cast<int>(br.sets.size()));
                }
            }
            force_aggressive_branching = false;
        }

        if (!used_balas) {
            int branch_var = -1;

            // Priority: branch on unfixed linking columns first
            if (decomp_state.enabled) {
                const auto &eff_dec = augmented_decisions.empty()
                                          ? branch_node.decisions
                                          : augmented_decisions;
                branch_var = pick_linking_branch_var(
                    decomp_state, eff_dec, fractional);
            }

            if (branch_var < 0) {
                if (use_reliability) {
                    int sb_lps = 0;
                    if (!adj_valid) {
                        cached_adj = build_adjacency(base);
                        adj_valid = true;
                    }
                    branch_var = reliability_branch_select(
                        pc_state, lp, base, branch_node, sol, fractional,
                        best_obj, config.reliability_eta, config.reliability_max_sb,
                        integ_tol, verbosity, sb_lps, &cached_adj);
                    total_lp_solves += sb_lps;
                } else {
                    branch_var = selector->select(sol.col_value, base.obj, fractional);
                }
            }
            if (branch_var < 0) continue;

            BranchNodeState child_zero;
            if (append_decision_if_consistent(branch_node, branch_var, 0, &child_zero)) {
                for (const auto &rc : rc_fixings)
                    child_zero.decisions.push_back(rc);
                if (!is_node_provably_infeasible(child_zero, base)) {
                    child_zero.parent_dual_bound = node_dual_bound;
                    child_zero.parent_dual_bound_raw = node_dual_bound_raw;
                    const int child_id = static_cast<int>(nodes.size());
                    nodes.push_back(std::move(child_zero));
                    node_bases[child_id] = node_basis;
                    frontier.push(child_id);
                }
            }

            BranchNodeState child_one;
            if (append_decision_if_consistent(branch_node, branch_var, 1, &child_one)) {
                for (const auto &rc : rc_fixings)
                    child_one.decisions.push_back(rc);
                if (!is_node_provably_infeasible(child_one, base)) {
                    child_one.parent_dual_bound = node_dual_bound;
                    child_one.parent_dual_bound_raw = node_dual_bound_raw;
                    const int child_id = static_cast<int>(nodes.size());
                    nodes.push_back(std::move(child_one));
                    node_bases[child_id] = node_basis;
                    frontier.push(child_id);
                }
            }
        }

        // Update global dual bound from heap top: O(1) instead of O(|frontier|)
        if (!frontier.empty()) {
            global_dual_bound = frontier.top_bound();
            global_dual_bound_raw = frontier.min_raw_bound;
        }

        // Stagnation control
        if (gap_stagnation_window > 0 && std::isfinite(best_obj)) {
            const double current_gap = compute_mip_gap(best_obj, global_dual_bound);
            if (std::isfinite(current_gap) && current_gap < best_mip_gap_seen - 1e-8) {
                best_mip_gap_seen = current_gap;
                node_at_last_gap_improvement = processed_nodes;
            }

            if (processed_nodes - node_at_last_gap_improvement >= gap_stagnation_window) {
                node_at_last_gap_improvement = processed_nodes;

                // Dynamic decomposition activation after ~150s of BnB
                if (!decomp_state.enabled && decomp_dynamic_eligible) {
                    const double bnb_elapsed = std::chrono::duration<double>(
                        Clock::now() - start_time).count();
                    if (bnb_elapsed > 150.0) {
                        decomp_state = compute_decomposition(
                            base, config.decomposition_max_linkers, verbosity);
                        decomp_dynamic_eligible = false;
                    }
                }

                // Track frontier trend with hysteresis:
                // Enter proving phase after 5 consecutive shrinks (sustained trend).
                // Exit only after 3 consecutive growths (avoid oscillation).
                if (frontier_at_last_stagnation > 0) {
                    if (frontier.size() < frontier_at_last_stagnation) {
                        frontier_shrink_streak = std::max(1, frontier_shrink_streak + 1);
                    } else {
                        frontier_shrink_streak = std::min(-1, frontier_shrink_streak - 1);
                    }
                }
                if (!proving_phase && frontier_shrink_streak >= 5) {
                    proving_phase = true;
                    if (verbosity >= 2) {
                        const double t = std::chrono::duration<double>(
                            Clock::now() - start_time).count();
                        fprintf(stderr,
                            "  [%8.3fs] Proving phase: frontier shrinking, "
                            "disabling expensive heuristics\n", t);
                    }
                }
                if (proving_phase && frontier_shrink_streak <= -3) {
                    proving_phase = false;
                    if (verbosity >= 2) {
                        const double t = std::chrono::duration<double>(
                            Clock::now() - start_time).count();
                        fprintf(stderr,
                            "  [%8.3fs] Exiting proving phase: frontier growing\n", t);
                    }
                }
                const bool frontier_shrinking = (frontier_shrink_streak >= 2);
                frontier_at_last_stagnation = frontier.size();

                // Mid-BnB cuts (skip if frontier is naturally shrinking)
                cut_accumulator += config.mid_bnb_cut_frequency;
                if (cut_accumulator >= 1.0 && !frontier_shrinking &&
                    config.cuts_enabled && sol.optimal) {
                    cut_accumulator -= 1.0;

                    // Save pre-cut state for rollback
                    const double pre_cut_dual = global_dual_bound;
                    const int pre_cut_nrows = base.nrows;
                    const int pre_cut_nnz = base.nnz;
                    const auto pre_cut_csr_inds = base.csr_inds;
                    const auto pre_cut_csr_offs = base.csr_offs;
                    const auto pre_cut_csr_vals = base.csr_vals;
                    const auto pre_cut_rhs = base.rhs;
                    const auto pre_cut_base_cuts = base.base_cuts;

                    int total_cuts = 0;
                    for (int round = 0; round < config.mid_bnb_cut_rounds; ++round) {
                        std::vector<CutConstraint> round_cuts;
                        for (const auto &sep : cut_separators) {
                            auto sep_cuts = sep->separate(sol.col_value, sol.row_dual,
                                                          base, base.ncols, integ_tol, best_obj);
                            for (auto &c : sep_cuts) {
                                if (static_cast<int>(round_cuts.size()) >= config.max_cuts_per_round) break;
                                round_cuts.push_back(std::move(c));
                            }
                            if (static_cast<int>(round_cuts.size()) >= config.max_cuts_per_round) break;
                        }
                        if (round_cuts.empty()) break;
                        append_cuts_to_base_model(base, round_cuts);
                        total_cuts += static_cast<int>(round_cuts.size());
                        // Use incremental add_cuts to preserve LP basis
                        lp.restore_base_state();
                        lp.add_cuts(round_cuts);
                        lp.apply_decisions(branch_node.decisions);
                        LpSolution cut_sol = lp.solve();
                        ++total_lp_solves;
                        if (!cut_sol.solved || !cut_sol.optimal) break;
                        lp.save_basis();
                        sol = cut_sol;
                    }

                    if (total_cuts > 0) {
                        // Rebuild to sync base model with LP
                        lp.rebuild_model(base);
                        adj_valid = false; orbital_valid = false;

                        // Validate: check if cuts improved the tightened dual
                        double post_cut_dual = sol.optimal ? sol.dual_obj : pre_cut_dual;
                        if (obj_is_integral && std::isfinite(post_cut_dual))
                            post_cut_dual = tighten_dual_bound(post_cut_dual, integ_tol);
                        if (post_cut_dual <= pre_cut_dual + tol) {
                            // Cuts didn't help — rollback
                            base.nrows = pre_cut_nrows;
                            base.nnz = pre_cut_nnz;
                            base.csr_inds = pre_cut_csr_inds;
                            base.csr_offs = pre_cut_csr_offs;
                            base.csr_vals = pre_cut_csr_vals;
                            base.rhs = pre_cut_rhs;
                            base.base_cuts = pre_cut_base_cuts;
                            lp.rebuild_model(base);
                            if (verbosity >= 3) {
                                const double t = std::chrono::duration<double>(Clock::now() - start_time).count();
                                fprintf(stderr, "  [%8.3fs] Stagnation: %d mid-BnB cuts undone (no dual improvement)\n",
                                        t, total_cuts);
                            }
                        } else {
                            if (verbosity >= 3) {
                                const double t = std::chrono::duration<double>(Clock::now() - start_time).count();
                                fprintf(stderr, "  [%8.3fs] Stagnation: added %d mid-BnB cuts (dual %.8f -> %.8f)\n",
                                        t, total_cuts, pre_cut_dual, post_cut_dual);
                            }
                        }
                    }
                }

                // Aggressive Balas branching (independent of cuts)
                balas_accumulator += config.aggressive_balas_frequency;
                if (balas_accumulator >= 1.0 && config.balas_enabled) {
                    balas_accumulator -= 1.0;
                    force_aggressive_branching = true;
                    if (verbosity >= 3) {
                        const double t = std::chrono::duration<double>(Clock::now() - start_time).count();
                        fprintf(stderr, "  [%8.3fs] Stagnation: aggressive Balas branching\n", t);
                    }
                }
            }
        }
    }

    // ================================================================
    // Final reporting
    // ================================================================
    // Recompute global dual bound
    {
        double new_bound = std::numeric_limits<double>::infinity();
        double new_bound_raw = std::numeric_limits<double>::infinity();
        for (const int idx : frontier.heap) {
            new_bound = std::min(new_bound, nodes[static_cast<size_t>(idx)].parent_dual_bound);
            new_bound_raw = std::min(new_bound_raw, nodes[static_cast<size_t>(idx)].parent_dual_bound_raw);
        }
        if (std::isfinite(new_bound)) global_dual_bound = new_bound;
        else if (frontier.empty() && std::isfinite(best_obj)) global_dual_bound = best_obj;
        if (std::isfinite(new_bound_raw)) global_dual_bound_raw = new_bound_raw;
        else if (frontier.empty() && std::isfinite(best_obj)) global_dual_bound_raw = best_obj;
    }

    const auto end_time = Clock::now();
    result.wall_time = std::chrono::duration<double>(end_time - start_time).count();
    result.nodes_processed = processed_nodes;
    result.lp_solves = total_lp_solves;

    // Add back the cost of columns fixed during preprocessing
    best_obj += fixed_preprocess_cost;
    global_dual_bound += fixed_preprocess_cost;
    global_dual_bound_raw += fixed_preprocess_cost;

    // Include fixed columns in the solution
    if (!fixed_original_cols.empty()) {
        if (best_solution.empty())
            best_solution.assign(static_cast<size_t>(ncols_input), 0.0);
        for (int c : fixed_original_cols) {
            if (c >= 0 && c < ncols_input)
                best_solution[static_cast<size_t>(c)] = 1.0;
        }
    }

    if (std::isfinite(best_obj)) {
        result.primal_obj = best_obj;
        result.solution = best_solution;

        if ((frontier_exhausted || gap_tolerance_reached) &&
            !hard_time_limit_reached &&
            processed_nodes < config.max_nodes) {
            result.dual_obj = best_obj;
            result.mip_gap = 0.0;
            result.status = "Optimal";
        } else {
            result.dual_obj = global_dual_bound;
            result.mip_gap = compute_mip_gap(best_obj, global_dual_bound);
            result.status = hard_time_limit_reached ? "TimeLimit" : "NodeLimit";
        }
    } else {
        result.primal_obj = std::numeric_limits<double>::infinity();
        result.dual_obj = global_dual_bound;
        result.mip_gap = std::numeric_limits<double>::infinity();
        result.status = "NoSolution";
    }

    return result;
}

} // namespace scpsol
