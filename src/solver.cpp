#include "solver.h"
#include "balas.h"
#include "bnb.h"
#include "compact.h"
#include "cuts.h"
#include "decomposition.h"
#include "diving.h"
#include "heuristics.h"
#include "lagrangian.h"
#include "lp.h"
#include "preprocessor.h"
#include "reliability.h"
#include "rins.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>

namespace scpsol {
namespace {

using Clock = std::chrono::steady_clock;
constexpr double kInf = std::numeric_limits<double>::infinity();

// ======================================================================
// Symmetry breaking: orbital fixing
// ======================================================================
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

OrbitalFixingState build_orbital_state(const BaseRelaxationModel &base) {
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
std::vector<BranchDecision> compute_orbital_fixings(
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
                    std::vector<int> cand_residual;
                    for (int r : adj.rows_by_col[static_cast<size_t>(cand)])
                        if (r >= 0 && r < nrows && !row_decided[static_cast<size_t>(r)])
                            cand_residual.push_back(r);
                    if (cand_residual == ref_residual) {
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

// ======================================================================
// Best-bound frontier: nodes ordered by ascending parent_dual_bound.
// Uses a vector managed as a min-heap for efficient best-bound selection.
// ======================================================================
struct BestBoundFrontier {
    std::vector<int> heap;
    const std::vector<BranchNodeState> *nodes_ptr = nullptr;
    double min_raw_bound = kInf;

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

    // Drop every node whose bound cannot beat `cutoff` (reduced-space value).
    void prune(double cutoff, double tol, int verbosity,
               std::unordered_map<int, HighsBasis> *bases = nullptr) {
        const double prune_bound = cutoff - tol;
        size_t before = heap.size();
        std::vector<int> surviving;
        surviving.reserve(heap.size());
        min_raw_bound = kInf;
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
        min_raw_bound = kInf;
        for (const int idx : heap) {
            const double raw = (*nodes_ptr)[static_cast<size_t>(idx)].parent_dual_bound_raw;
            if (raw < min_raw_bound) min_raw_bound = raw;
        }
    }
};

// Apply a column renumbering produced by a mid-BnB reduction to every open
// node, the pseudocosts and (optionally) the node currently being processed.
void apply_mid_bnb_reduction(
    const ModelReductionResult &reduction,
    BaseRelaxationModel &base,
    BestBoundFrontier &frontier,
    std::vector<BranchNodeState> &nodes,
    PseudocostState *pc_state,
    BranchNodeState *current_node,
    bool *current_alive) {
    base.build_transpose();
    std::vector<int> surviving;
    for (const int idx : frontier.heap) {
        if (remap_branch_node(nodes[static_cast<size_t>(idx)], reduction.old_to_new))
            surviving.push_back(idx);
    }
    frontier.heap = std::move(surviving);
    frontier.rebuild_heap();
    frontier.recompute_min_raw();
    if (pc_state) pc_state->remap(reduction.old_to_new, base.ncols);
    if (current_node && current_alive && *current_alive)
        *current_alive = remap_branch_node(*current_node, reduction.old_to_new);
}

int mid_bnb_column_removal(
    BaseRelaxationModel &base,
    double cutoff,
    double tol,
    BestBoundFrontier &frontier,
    std::vector<BranchNodeState> &nodes,
    int verbosity,
    PseudocostState *pc_state,
    BranchNodeState *current_node,
    bool *current_alive) {
    ModelReductionResult reduction = reduce_base_model(base, cutoff, tol);
    if (reduction.columns_removed <= 0) return 0;
    if (verbosity >= 3)
        fprintf(stderr, "           Mid-BnB reduction: %d cols removed, %d remaining\n",
                reduction.columns_removed, base.ncols);
    apply_mid_bnb_reduction(reduction, base, frontier, nodes, pc_state, current_node, current_alive);
    return reduction.columns_removed;
}

int mid_bnb_budget_pruning(
    BaseRelaxationModel &base,
    double cutoff,
    double tol,
    double preprocess_time_limit,
    BestBoundFrontier &frontier,
    std::vector<BranchNodeState> &nodes,
    int verbosity,
    PseudocostState *pc_state,
    BranchNodeState *current_node,
    bool *current_alive) {
    ModelReductionResult reduction = reduce_base_model_budget_pruning(
        base, cutoff, tol, preprocess_time_limit);
    if (reduction.columns_removed <= 0) return 0;
    if (verbosity >= 3)
        fprintf(stderr, "           Mid-BnB budget pruning: %d cols removed, %d remaining\n",
                reduction.columns_removed, base.ncols);
    apply_mid_bnb_reduction(reduction, base, frontier, nodes, pc_state, current_node, current_alive);
    return reduction.columns_removed;
}

// ======================================================================
// Incumbent and shared solver state
// ======================================================================

// The incumbent always lives in the column space of the input instance and
// includes every column fixed during preprocessing. It is never "invalidated"
// by model reductions: reductions only remove columns that cannot appear in
// a strictly better solution, which the stored incumbent does not need.
struct Incumbent {
    double obj = kInf;              // total cost
    std::vector<double> solution;   // one 0/1 entry per input column
    std::string source = "none";
};

struct SolverRun {
    const SolverConfig &config;
    const int verbosity;
    const double tol;
    const double integ_tol;
    const Clock::time_point start_time;
    const int ncols_input;
    const bool obj_is_integral;

    // Reduced instance manipulated by root preprocessing.
    WorkingModel wm;

    Incumbent inc;
    // Cost of the columns permanently fixed to 1 during preprocessing. The
    // reduced models never see these columns, so every bound they produce is
    // offset by this amount.
    double fixed_cost = 0.0;
    std::vector<int> fixed_original_cols;

    // Dual bounds in reduced space (add fixed_cost for the original space).
    double global_dual_bound = -kInf;     // tightened
    double global_dual_bound_raw = -kInf; // raw LP value

    BaseRelaxationModel base;
    LpSolver lp;
    LpSolution root_sol;
    int total_lp_solves = 0;
    std::vector<std::unique_ptr<IIntegerHeuristic>> heuristics;

    SolverRun(const ScpInstance &instance, const SolverConfig &cfg)
        : config(cfg),
          verbosity(cfg.verbosity),
          tol(cfg.feasibility_tol),
          integ_tol(cfg.integrality_tol),
          start_time(Clock::now()),
          ncols_input(instance.ncols),
          obj_is_integral(has_integer_objective(instance.costs, instance.ncols, cfg.integrality_tol)) {
        wm.nrows = instance.nrows;
        wm.ncols = instance.ncols;
        wm.csr_inds = instance.csr_indices;
        wm.csr_offs = instance.csr_offsets;
        wm.csr_vals = instance.csr_values;
        wm.obj = instance.costs;
        wm.active_to_input.resize(static_cast<size_t>(instance.ncols));
        std::iota(wm.active_to_input.begin(), wm.active_to_input.end(), 0);
        heuristics = make_integer_heuristics(cfg.heuristic_config);
    }

    double elapsed() const {
        return std::chrono::duration<double>(Clock::now() - start_time).count();
    }

    // Objective value (in reduced space) that a solution of the current
    // reduced model must beat to improve on the incumbent.
    double cutoff() const { return inc.obj - fixed_cost; }

    // Record columns permanently fixed to 1 by preprocessing.
    void register_fixed(double cost, const std::vector<int> &original_cols) {
        if (original_cols.empty() && cost == 0.0) return;
        fixed_cost += cost;
        fixed_original_cols.insert(fixed_original_cols.end(),
                                   original_cols.begin(), original_cols.end());
        if (std::isfinite(global_dual_bound)) global_dual_bound -= cost;
        if (std::isfinite(global_dual_bound_raw)) global_dual_bound_raw -= cost;
    }

    // Rebuild the LP base model from the working model.
    void build_base() {
        base.nrows = wm.nrows;
        base.ncols = wm.ncols;
        base.ncols_input = ncols_input;
        base.nnz = static_cast<int>(wm.csr_vals.size());
        base.csr_inds = wm.csr_inds;
        base.csr_offs = wm.csr_offs;
        base.csr_vals = wm.csr_vals;
        base.obj = wm.obj;
        base.rhs.assign(static_cast<size_t>(wm.nrows), 1.0);
        base.active_to_original = wm.active_to_input;
        base.base_cuts.clear();
        base.build_transpose();
    }

    // Adopt a solution of the current base model (reduced space) if it
    // improves on the incumbent. Returns true on improvement.
    bool try_adopt(double reduced_obj, const std::vector<double> &reduced_sol,
                   const std::string &source) {
        const double total = reduced_obj + fixed_cost;
        if (!(total < inc.obj - tol)) return false;

        inc.obj = total;
        inc.source = source;
        inc.solution.assign(static_cast<size_t>(ncols_input), 0.0);
        const int copy_cols = std::min(base.ncols, static_cast<int>(reduced_sol.size()));
        for (int j = 0; j < copy_cols; ++j) {
            if (reduced_sol[static_cast<size_t>(j)] > 0.5) {
                const int orig = base.active_to_original[static_cast<size_t>(j)];
                if (orig >= 0 && orig < ncols_input)
                    inc.solution[static_cast<size_t>(orig)] = 1.0;
            }
        }
        for (int c : fixed_original_cols) {
            if (c >= 0 && c < ncols_input)
                inc.solution[static_cast<size_t>(c)] = 1.0;
        }
        if (verbosity >= 2)
            fprintf(stderr, "* [%8.3fs] New incumbent: %.8f (from %s)\n",
                    elapsed(), inc.obj, source.c_str());
        return true;
    }

    // Run the primal heuristics on an LP solution of the base model.
    // Returns true if the incumbent improved.
    bool run_heuristics(const LpSolution &sol, const BranchNodeState &node,
                        const std::string &prefix, bool stop_at_first) {
        bool improved = false;
        for (const auto &h : heuristics) {
            auto hr = h->tryBuild(sol.col_value, sol.row_dual, base, node, integ_tol);
            if (hr.feasible && try_adopt(hr.objective, hr.solution, prefix + hr.name)) {
                improved = true;
                if (stop_at_first) break;
            }
        }
        return improved;
    }
};

// ======================================================================
// Root reductions on the working model
// ======================================================================

void apply_cost_reduction(SolverRun &R) {
    WorkingModel &wm = R.wm;
    const double cutoff = R.cutoff();
    if (!std::isfinite(cutoff) || wm.ncols <= 0) return;

    std::vector<char> keep(static_cast<size_t>(wm.ncols), 1);
    int removed = 0;
    for (int j = 0; j < wm.ncols; ++j) {
        if (wm.obj[static_cast<size_t>(j)] + R.tol >= cutoff) {
            keep[static_cast<size_t>(j)] = 0;
            ++removed;
        }
    }
    if (removed == 0) return;
    wm.remove_columns(keep);
    if (R.verbosity >= 3)
        fprintf(stderr, "  Cost reduction: %d cols removed, %d remaining\n", removed, wm.ncols);
}

void apply_budget_pruning(SolverRun &R) {
    WorkingModel &wm = R.wm;
    const double cutoff = R.cutoff();
    if (!std::isfinite(cutoff) || wm.ncols <= 0) return;

    ColumnPreprocessContext ctx = make_column_context(
        wm.nrows, wm.ncols, wm.csr_inds, wm.csr_offs, wm.csr_vals, wm.obj,
        R.tol, R.config.preprocess_time_limit, cutoff);
    auto rules = make_preprocess_rules("incumbent_budget");
    int total_removed = 0;
    for (const auto &rule : rules) total_removed += rule->apply(ctx, R.tol);
    if (total_removed <= 0) return;

    wm.remove_columns(ctx.active);
    if (R.verbosity >= 3)
        fprintf(stderr, "  Budget pruning: %d cols removed, %d remaining\n", total_removed, wm.ncols);
}

void apply_dominance_reduction(SolverRun &R, const std::string &rules_config) {
    WorkingModel &wm = R.wm;
    if (wm.ncols <= 0) return;

    ColumnPreprocessContext ctx = make_column_context(
        wm.nrows, wm.ncols, wm.csr_inds, wm.csr_offs, wm.csr_vals, wm.obj,
        R.tol, R.config.preprocess_time_limit, kInf);
    auto rules = make_preprocess_rules(rules_config);
    int total_removed = 0;
    for (const auto &rule : rules) total_removed += rule->apply(ctx, R.tol);
    if (total_removed <= 0) return;

    const int before = wm.ncols;
    wm.remove_columns(ctx.active);
    if (R.verbosity >= 3)
        fprintf(stderr, "  Dominance reduction: %d -> %d cols\n", before, wm.ncols);
}

void apply_dominance_finder(SolverRun &R) {
    WorkingModel &wm = R.wm;
    if (wm.ncols <= 0 || wm.nrows <= 0) return;

    ColumnPreprocessContext ctx = make_column_context(
        wm.nrows, wm.ncols, wm.csr_inds, wm.csr_offs, wm.csr_vals, wm.obj,
        R.tol, R.config.preprocess_time_limit, kInf);
    const int removed = dominance_finder(ctx, R.tol, 0.5, R.verbosity);
    if (removed <= 0) return;

    wm.remove_columns(ctx.active);
    if (R.verbosity >= 2)
        fprintf(stderr, "  Dominance Finder: %d cols removed, %d remaining\n", removed, wm.ncols);
}

// Essential-column fixing, row domination and probing. Returns true if the
// model changed.
bool apply_row_reduction(SolverRun &R) {
    WorkingModel &wm = R.wm;
    auto rr = row_reduce(wm.nrows, wm.ncols, wm.csr_inds, wm.csr_offs, wm.csr_vals,
                         wm.obj, wm.active_to_input,
                         R.tol, R.config.preprocess_time_limit, R.verbosity);
    const bool changed = (rr.cols_fixed > 0 || rr.rows_removed > 0);
    if (changed) R.register_fixed(rr.fixed_cost, rr.fixed_original_cols);
    return changed;
}

// Re-run the cheap reductions until nothing changes (bounded number of rounds).
void iterate_reductions(SolverRun &R, const char *tag) {
    const int max_extra_rounds = 9;
    for (int round = 0; round < max_extra_rounds; ++round) {
        const int ncols_start = R.wm.ncols;
        const int nrows_start = R.wm.nrows;
        apply_cost_reduction(R);
        apply_dominance_reduction(R, R.config.preprocess_rules);
        apply_row_reduction(R);
        if (R.wm.ncols == ncols_start && R.wm.nrows == nrows_start) break;
        if (R.verbosity >= 3)
            fprintf(stderr, "  %s round %d: %d/%d -> %d/%d (cols/rows)\n",
                    tag, round + 2, ncols_start, nrows_start, R.wm.ncols, R.wm.nrows);
    }
}

// ======================================================================
// Phase 1: greedy heuristic
// ======================================================================
void run_greedy(SolverRun &R) {
    if (R.verbosity >= 2) fprintf(stderr, "Phase 1: Greedy set cover heuristic\n");
    const WorkingModel &wm = R.wm;
    auto greedy = greedy_set_cover_heuristic(wm.nrows, wm.ncols, wm.csr_inds, wm.csr_offs,
                                             wm.csr_vals, wm.obj.data());
    if (!greedy.feasible) return;

    R.inc.obj = greedy.objective;
    R.inc.source = "greedy";
    R.inc.solution.assign(static_cast<size_t>(R.ncols_input), 0.0);
    for (int col : greedy.selected_columns) {
        const int input_col = wm.active_to_input[static_cast<size_t>(col)];
        if (input_col >= 0 && input_col < R.ncols_input)
            R.inc.solution[static_cast<size_t>(input_col)] = 1.0;
    }
    if (R.verbosity >= 2)
        fprintf(stderr, "  Greedy incumbent: %.12g\n", R.inc.obj);
}

// ======================================================================
// Phase 2: cost + budget + dominance + row reduction (iterated)
// ======================================================================
void run_root_reductions(SolverRun &R) {
    const int cols_before = R.wm.ncols;

    // Step 1: cheap reductions (cost + budget)
    apply_cost_reduction(R);
    apply_budget_pruning(R);

    // Step 2: greedy multi-column dominance (fast, removes bulk of dominated cols)
    apply_dominance_finder(R);

    // Step 3: exact pairwise/triplet dominance on the reduced set
    apply_dominance_reduction(R, "cost_driven");
    apply_dominance_reduction(R, R.config.preprocess_rules);

    // Step 4: row reduction (essential columns, row domination, probing)
    if (apply_row_reduction(R))
        iterate_reductions(R, "Preprocess");

    if (R.wm.ncols < cols_before && R.verbosity >= 2)
        fprintf(stderr, "  Pre-LP reduction: cols %d -> %d, rows %d\n",
                cols_before, R.wm.ncols, R.wm.nrows);
}

// Result when preprocessing alone covered every row with essential columns.
SolverResult trivial_result(const SolverRun &R) {
    SolverResult result;
    result.wall_time = R.elapsed();
    result.primal_obj = R.fixed_cost;
    result.dual_obj = R.fixed_cost;
    result.mip_gap = 0.0;
    result.nodes_processed = 0;
    result.lp_solves = 0;
    result.status = "Optimal";
    result.solution.assign(static_cast<size_t>(R.ncols_input), 0.0);
    for (int c : R.fixed_original_cols) {
        if (c >= 0 && c < R.ncols_input)
            result.solution[static_cast<size_t>(c)] = 1.0;
    }
    if (R.verbosity >= 1)
        fprintf(stderr, "All rows covered by essential columns (cost %.12g)\n", R.fixed_cost);
    return result;
}

// ======================================================================
// Phase 3: root LP relaxation and Lagrangian heuristic
// ======================================================================
void solve_root_lp(SolverRun &R) {
    if (R.verbosity >= 2) fprintf(stderr, "Phase 3: Root LP relaxation\n");

    R.build_base();
    R.lp.build_model(R.base);

    // Root LP: use IPM for large models, simplex otherwise
    const bool use_ipm_root = (R.config.root_lp_solver == "ipm") ||
        (R.config.root_lp_solver == "auto" && R.base.ncols >= 500 && R.base.nrows >= 200);
    if (use_ipm_root) {
        if (R.verbosity >= 2)
            fprintf(stderr, "  Root LP solver: IPM (barrier) with crossover\n");
        R.root_sol = R.lp.solve_ipm();
    } else {
        R.root_sol = R.lp.solve();
    }
    ++R.total_lp_solves;

    if (!R.root_sol.solved || !R.root_sol.optimal) {
        if (R.verbosity >= 1)
            fprintf(stderr, "  Root LP did not converge\n");
        return;
    }
    R.lp.save_basis();

    // Heuristics on the root LP
    BranchNodeState root_node;
    R.run_heuristics(R.root_sol, root_node, "root_", false);

    // Integral root LP
    if (static_cast<int>(R.root_sol.col_value.size()) >= R.base.ncols &&
        is_binary_integral(R.root_sol.col_value, R.base.ncols, R.integ_tol))
        R.try_adopt(R.root_sol.primal_obj, R.root_sol.col_value, "root_lp_exact");

    double root_dual = R.root_sol.dual_obj;
    R.global_dual_bound_raw = root_dual;
    if (R.obj_is_integral)
        root_dual = tighten_dual_bound(root_dual, R.integ_tol);
    R.global_dual_bound = root_dual;

    if (R.verbosity >= 2)
        fprintf(stderr, "  Root LP: primal=%.12g dual=%.12g\n",
                R.root_sol.primal_obj, R.root_sol.dual_obj);

    // Phase 3.5: Lagrangian relaxation (subgradient optimization) as a heuristic
    const int lagr_max_iters = 500;
    const double lagr_time_limit = 5.0;
    auto lagr = lagrangian_relaxation(R.base, R.cutoff(), R.root_sol.row_dual,
                                      lagr_max_iters, lagr_time_limit, R.verbosity);
    if (!lagr.best_solution.empty())
        R.try_adopt(lagr.best_heuristic_obj, lagr.best_solution, "lagrangian");
}

// ======================================================================
// Phase 4: iterative reduced-cost fixing (fix-resolve-fix cycles)
// ======================================================================
void run_rc_fixing(SolverRun &R) {
    if (!R.root_sol.solved || !R.root_sol.optimal || !std::isfinite(R.cutoff())) return;

    WorkingModel &wm = R.wm;
    const int max_rc_iters = 10;
    int total_rc_removed = 0;
    for (int rc_iter = 0; rc_iter < max_rc_iters; ++rc_iter) {
        const double gap = R.cutoff() - R.root_sol.dual_obj;
        if (gap <= R.tol) break;

        const int n = std::min(wm.ncols, static_cast<int>(R.root_sol.col_value.size()));
        auto rcosts = compute_reduced_costs(wm.obj, R.root_sol.row_dual, R.base, n);
        std::vector<char> keep(static_cast<size_t>(wm.ncols), 1);
        int rc_removed = 0;
        for (int j = 0; j < n; ++j) {
            if (R.root_sol.col_value[static_cast<size_t>(j)] < R.tol &&
                rcosts[static_cast<size_t>(j)] > gap + R.tol) {
                keep[static_cast<size_t>(j)] = 0;
                ++rc_removed;
            }
        }
        if (rc_removed == 0) break;

        wm.remove_columns(keep);
        total_rc_removed += rc_removed;
        if (R.verbosity >= 3)
            fprintf(stderr, "  RC fixing iter %d: %d cols removed, %d remaining\n",
                    rc_iter + 1, rc_removed, wm.ncols);

        // Rebuild base model and re-solve LP
        R.build_base();
        R.lp.rebuild_model(R.base);
        R.root_sol = R.lp.solve();
        ++R.total_lp_solves;
        if (!R.root_sol.solved || !R.root_sol.optimal) break;
        R.lp.save_basis();

        // Heuristics on the tighter LP (may improve the incumbent and the gap)
        BranchNodeState rc_node;
        R.run_heuristics(R.root_sol, rc_node, "rcfix_", false);
    }
    if (total_rc_removed > 0 && R.verbosity >= 2)
        fprintf(stderr, "  Iterative RC fixing: %d cols removed total\n", total_rc_removed);
}

// ======================================================================
// Phase 5: LP-based probing (fix x_j = 1 when x_j = 0 cannot improve)
// ======================================================================
void run_lp_probing(SolverRun &R) {
    if (!R.root_sol.solved || !R.root_sol.optimal || !std::isfinite(R.cutoff()) ||
        R.wm.ncols > 800)
        return;

    WorkingModel &wm = R.wm;
    const auto probe_start = Clock::now();
    const double probe_time_limit = 10.0;

    // Only probe columns with x_j > 0 in the LP (forcing a zero to 0 is a no-op),
    // most impactful first.
    std::vector<int> probe_order;
    for (int j = 0; j < wm.ncols; ++j) {
        if (R.root_sol.col_value[static_cast<size_t>(j)] > R.tol)
            probe_order.push_back(j);
    }
    std::sort(probe_order.begin(), probe_order.end(), [&](int a, int b) {
        return R.root_sol.col_value[static_cast<size_t>(a)] >
               R.root_sol.col_value[static_cast<size_t>(b)];
    });

    std::vector<char> probe_fix_one(static_cast<size_t>(wm.ncols), 0);
    int probe_fixed = 0;
    HighsBasis probe_basis = R.lp.get_basis();

    for (int j : probe_order) {
        const double elapsed = std::chrono::duration<double>(Clock::now() - probe_start).count();
        if (elapsed > probe_time_limit) break;

        R.lp.restore_base_state();
        R.lp.set_basis(probe_basis);
        BranchDecision fix_zero{j, 0};
        R.lp.apply_decisions({fix_zero});
        auto probe_sol = R.lp.solve();
        ++R.total_lp_solves;

        // Only a proven result may fix a column: an LP that failed to solve
        // says nothing about x_j = 0.
        if (probe_sol.infeasible ||
            (probe_sol.optimal && probe_sol.primal_obj > R.cutoff() + R.tol)) {
            probe_fix_one[static_cast<size_t>(j)] = 1;
            ++probe_fixed;
        }
    }

    // Restore LP state
    R.lp.restore_base_state();
    R.lp.set_basis(probe_basis);

    if (probe_fixed == 0) return;
    if (R.verbosity >= 2)
        fprintf(stderr, "  LP probing: %d cols fixed to 1\n", probe_fixed);

    double probe_fixed_cost = 0.0;
    std::vector<int> probe_fixed_cols;
    std::vector<char> keep_row(static_cast<size_t>(wm.nrows), 1);
    std::vector<char> keep_col(static_cast<size_t>(wm.ncols), 1);
    for (int j = 0; j < wm.ncols; ++j) {
        if (!probe_fix_one[static_cast<size_t>(j)]) continue;
        keep_col[static_cast<size_t>(j)] = 0;
        probe_fixed_cost += wm.obj[static_cast<size_t>(j)];
        probe_fixed_cols.push_back(wm.active_to_input[static_cast<size_t>(j)]);
        for (const auto &entry : R.base.cols_to_rows[static_cast<size_t>(j)]) {
            if (entry.row >= 0 && entry.row < wm.nrows && entry.val > 0.0)
                keep_row[static_cast<size_t>(entry.row)] = 0;
        }
    }

    // Drop the fixed columns and the rows they cover
    wm.remove_columns(keep_col);
    wm.remove_rows(keep_row);
    R.register_fixed(probe_fixed_cost, probe_fixed_cols);

    // Rebuild and re-solve for heuristics on the reduced model
    R.build_base();
    R.lp.rebuild_model(R.base);
    R.root_sol = R.lp.solve();
    ++R.total_lp_solves;
    if (R.root_sol.solved && R.root_sol.optimal) {
        R.lp.save_basis();
        BranchNodeState probe_node;
        R.run_heuristics(R.root_sol, probe_node, "probe_", false);
    }
}

// ======================================================================
// Phases 4-5: everything between the root LP and the BnB base model
// ======================================================================
void run_post_lp_phase(SolverRun &R) {
    const int cols_before = R.wm.ncols;

    run_rc_fixing(R);
    run_lp_probing(R);

    // Post-LP reductions with the (possibly improved) incumbent
    apply_cost_reduction(R);
    apply_budget_pruning(R);
    apply_dominance_reduction(R, R.config.preprocess_rules);
    if (apply_row_reduction(R))
        iterate_reductions(R, "Post-LP");

    if (R.wm.ncols < cols_before && R.verbosity >= 2)
        fprintf(stderr, "  Post-LP reduction: cols %d -> %d, rows %d\n",
                cols_before, R.wm.ncols, R.wm.nrows);
}

// ======================================================================
// Phase 6: BnB base model, root cuts, root Balas cuts, post-cut pruning
// ======================================================================
void run_root_cut_rounds(SolverRun &R) {
    if (!R.config.cuts_enabled || R.config.cut_rounds_root <= 0 ||
        !R.root_sol.solved || !R.root_sol.optimal)
        return;

    BaseRelaxationModel &base = R.base;
    auto separators = make_cut_separators();
    LpSolution cut_sol = R.root_sol;
    int root_cuts_added = 0;

    // Save pre-cut model state to restore if cuts don't help
    const double pre_cut_tightened_dual = R.global_dual_bound;
    const int pre_cut_nrows = base.nrows;
    const int pre_cut_nnz = base.nnz;
    const auto pre_cut_csr_inds = base.csr_inds;
    const auto pre_cut_csr_offs = base.csr_offs;
    const auto pre_cut_csr_vals = base.csr_vals;
    const auto pre_cut_rhs = base.rhs;
    const auto pre_cut_base_cuts = base.base_cuts;

    for (int round = 0; round < R.config.cut_rounds_root; ++round) {
        // Integral LP after cuts
        if (static_cast<int>(cut_sol.col_value.size()) >= base.ncols &&
            is_binary_integral(cut_sol.col_value, base.ncols, R.integ_tol) &&
            R.try_adopt(cut_sol.primal_obj, cut_sol.col_value, "cut_round_exact")) {
            if (R.verbosity >= 2)
                fprintf(stderr, "  Cut round %d: LP integral, incumbent %.12g\n", round + 1, R.inc.obj);
            break;
        }

        BranchNodeState empty_branch;
        R.run_heuristics(cut_sol, empty_branch, "cut_", false);

        // Update dual bound
        if (cut_sol.optimal && std::isfinite(cut_sol.dual_obj)) {
            double cut_dual = cut_sol.dual_obj;
            if (cut_dual > R.global_dual_bound_raw) R.global_dual_bound_raw = cut_dual;
            if (R.obj_is_integral) cut_dual = tighten_dual_bound(cut_dual, R.integ_tol);
            if (cut_dual > R.global_dual_bound) R.global_dual_bound = cut_dual;
        }

        // Separate
        std::vector<CutConstraint> round_cuts;
        for (const auto &sep : separators) {
            auto sep_cuts = sep->separate(cut_sol.col_value, cut_sol.row_dual,
                                          base, base.ncols, R.integ_tol, R.cutoff());
            for (auto &c : sep_cuts) {
                if (static_cast<int>(round_cuts.size()) >= R.config.max_cuts_per_round) break;
                round_cuts.push_back(std::move(c));
            }
            if (static_cast<int>(round_cuts.size()) >= R.config.max_cuts_per_round) break;
        }
        if (round_cuts.empty()) {
            if (R.verbosity >= 3)
                fprintf(stderr, "  Cut round %d: no violated cuts\n", round + 1);
            break;
        }

        append_cuts_to_base_model(base, round_cuts);
        root_cuts_added += static_cast<int>(round_cuts.size());
        if (R.verbosity >= 3)
            fprintf(stderr, "  Cut round %d: added %d cuts (model %d rows)\n",
                    round + 1, static_cast<int>(round_cuts.size()), base.nrows);

        // Rebuild LP with cuts
        R.lp.rebuild_model(base);
        cut_sol = R.lp.solve();
        ++R.total_lp_solves;
        if (!cut_sol.solved || !cut_sol.optimal) break;
        R.lp.save_basis();
    }

    // For integer objectives, cuts that don't push the LP past the next
    // integer provide no benefit and can disrupt branching topology.
    if (root_cuts_added > 0 && R.global_dual_bound <= pre_cut_tightened_dual + R.tol) {
        base.nrows = pre_cut_nrows;
        base.nnz = pre_cut_nnz;
        base.csr_inds = pre_cut_csr_inds;
        base.csr_offs = pre_cut_csr_offs;
        base.csr_vals = pre_cut_csr_vals;
        base.rhs = pre_cut_rhs;
        base.base_cuts = pre_cut_base_cuts;
        R.lp.rebuild_model(base);
        if (R.verbosity >= 2)
            fprintf(stderr, "  Root cuts undone: no tightened dual improvement\n");
    } else if (root_cuts_added > 0 && R.verbosity >= 2) {
        fprintf(stderr, "  Root cuts: %d total (dual %.12g -> %.12g)\n",
                root_cuts_added, pre_cut_tightened_dual, R.global_dual_bound);
    }
}

void run_root_balas_cuts(SolverRun &R) {
    if (!R.config.balas_enabled || !R.root_sol.solved || !R.root_sol.optimal) return;

    BaseRelaxationModel &base = R.base;
    LpSolution balas_sol = R.lp.solve();
    ++R.total_lp_solves;
    if (!balas_sol.solved || !balas_sol.optimal) return;

    auto rcosts = compute_reduced_costs(base.obj, balas_sol.row_dual, base, base.ncols);
    auto br = balas_branch_generate(balas_sol.col_value, balas_sol.row_dual, rcosts,
                                    R.cutoff(), base, base.ncols, R.config.balas_max_branches,
                                    R.integ_tol, false);
    if (br.sets.empty()) return;

    std::vector<CutConstraint> balas_cuts;
    for (const auto &R_k : br.sets) {
        // Only add cover cuts that are violated by the LP
        double lhs_val = 0.0;
        for (const int j : R_k) {
            if (j >= 0 && j < base.ncols)
                lhs_val += balas_sol.col_value[static_cast<size_t>(j)];
        }
        if (lhs_val >= 1.0 - R.integ_tol) continue;

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
    if (balas_cuts.empty()) return;

    append_cuts_to_base_model(base, balas_cuts);
    R.lp.rebuild_model(base);
    if (R.verbosity >= 2)
        fprintf(stderr, "  Root Balas cover cuts: %d added (from %d sets)\n",
                static_cast<int>(balas_cuts.size()), static_cast<int>(br.sets.size()));
}

void prepare_bnb_model(SolverRun &R) {
    R.build_base();
    R.lp.rebuild_model(R.base);

    if (R.root_sol.solved && R.root_sol.optimal) {
        // Re-solve to get a warm basis for the reduced model
        R.root_sol = R.lp.solve();
        ++R.total_lp_solves;
        if (R.root_sol.solved && R.root_sol.optimal) R.lp.save_basis();
    }

    if (R.verbosity >= 2) {
        const BaseRelaxationModel &base = R.base;
        const double density = (base.nrows > 0 && base.ncols > 0)
            ? 100.0 * base.nnz / (static_cast<double>(base.nrows) * base.ncols)
            : 0.0;
        fprintf(stderr, "BnB base model: %d rows x %d cols, nnz=%d (%.1f%% dense)\n",
                base.nrows, base.ncols, base.nnz, density);
    }

    run_root_cut_rounds(R);
    run_root_balas_cuts(R);

    // Post-cut budget pruning
    if (std::isfinite(R.cutoff())) {
        ModelReductionResult budget_red = reduce_base_model_budget_pruning(
            R.base, R.cutoff(), R.tol, R.config.preprocess_time_limit);
        if (budget_red.columns_removed > 0) {
            if (R.verbosity >= 3)
                fprintf(stderr, "  Post-cut budget pruning: %d cols removed\n",
                        budget_red.columns_removed);
            R.base.build_transpose();
            R.lp.rebuild_model(R.base);
        }
    }
}

// ======================================================================
// Branch-and-bound
// ======================================================================
struct BnbOutcome {
    int processed_nodes = 0;
    bool gap_tolerance_reached = false;
    bool hard_time_limit_reached = false;
    bool frontier_exhausted = false;
};

BnbOutcome run_branch_and_bound(SolverRun &R) {
    BnbOutcome out;
    const SolverConfig &config = R.config;
    const int verbosity = R.verbosity;
    const double tol = R.tol;
    const double integ_tol = R.integ_tol;
    BaseRelaxationModel &base = R.base;
    LpSolver &lp = R.lp;

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
        root_state.parent_dual_bound = std::isfinite(R.global_dual_bound)
                                           ? R.global_dual_bound : -kInf;
        root_state.parent_dual_bound_raw = std::isfinite(R.global_dual_bound_raw)
                                               ? R.global_dual_bound_raw : -kInf;
        nodes.push_back(root_state);
    }

    BestBoundFrontier frontier;
    frontier.init(&nodes);
    frontier.push(0);

    int processed_nodes = 0;

    const int gap_stagnation_window = config.gap_stagnation_window;
    double best_mip_gap_seen = kInf;
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

    // Everything that has to happen after the incumbent improved: prune the
    // frontier, drop columns that can no longer be in an improving solution,
    // and keep the LP, adjacency and symmetry caches in sync.
    //
    // If columns were removed while `current` (the node being processed) is
    // still needed, it is remapped and re-enqueued with the given bound, and
    // the caller must stop processing it (its LP solution and basis refer to
    // the old column space). Returns true in that case.
    auto on_new_incumbent = [&](BranchNodeState *current, double cur_bound,
                                double cur_bound_raw) -> bool {
        node_at_last_gap_improvement = processed_nodes;
        frontier.prune(R.cutoff(), tol, verbosity, &node_bases);

        bool current_alive = (current != nullptr);
        bool changed = false;
        if (mid_bnb_column_removal(base, R.cutoff(), tol, frontier, nodes, verbosity,
                                   &pc_state, current, &current_alive) > 0)
            changed = true;
        if (mid_bnb_budget_pruning(base, R.cutoff(), tol, config.preprocess_time_limit,
                                   frontier, nodes, verbosity,
                                   &pc_state, current, &current_alive) > 0)
            changed = true;
        if (!changed) return false;

        lp.rebuild_model_keep_basis(base);
        adj_valid = false;
        orbital_valid = false;

        if (current && current_alive && cur_bound < R.cutoff() - tol) {
            BranchNodeState requeue = *current;
            requeue.parent_dual_bound = cur_bound;
            requeue.parent_dual_bound_raw = cur_bound_raw;
            const int id = static_cast<int>(nodes.size());
            nodes.push_back(std::move(requeue));
            frontier.push(id);
            if (verbosity >= 3)
                fprintf(stderr, "           Current node re-enqueued after column renumbering\n");
        }
        return true;
    };

    if (verbosity >= 2)
        fprintf(stderr, "Branch-and-bound started (max_nodes=%d)\n", config.max_nodes);

    while (processed_nodes < config.max_nodes) {
        // Time limit check
        const auto now = Clock::now();
        const double elapsed = std::chrono::duration<double>(now - R.start_time).count();
        if (config.time_limit_seconds > 0.0 && elapsed >= config.time_limit_seconds) {
            out.hard_time_limit_reached = true;
            if (verbosity >= 1)
                fprintf(stderr, "  [%8.3fs] Time limit reached\n", elapsed);
            break;
        }

        // Gap check (original space)
        if (std::isfinite(R.inc.obj) && std::isfinite(R.global_dual_bound)) {
            const double gap = compute_mip_gap(R.inc.obj, R.global_dual_bound + R.fixed_cost);
            if (std::isfinite(gap) && gap <= config.mip_gap_tol) {
                out.gap_tolerance_reached = true;
                if (verbosity >= 1)
                    fprintf(stderr, "  [%8.3fs] MIP gap %.8f%% within tolerance; optimal\n",
                            elapsed, gap * 100.0);
                break;
            }
        }

        // Log
        if (config.log_interval_seconds > 0.0) {
            const double bnb_elapsed = std::chrono::duration<double>(now - bnb_start).count();
            if (bnb_elapsed >= next_log_sec) {
                if (verbosity >= 2) {
                    const double dual_raw = R.global_dual_bound_raw + R.fixed_cost;
                    const double gap = compute_mip_gap(R.inc.obj, dual_raw);
                    fprintf(stderr, "  [%8.3fs] nodes=%6d frontier=%6zu lp=%6d incumbent=%.8f dual=%.8f gap=",
                            elapsed, processed_nodes, frontier.size(), R.total_lp_solves,
                            R.inc.obj, dual_raw);
                    if (std::isfinite(gap)) fprintf(stderr, "%.8f%%\n", gap * 100.0);
                    else fprintf(stderr, "inf\n");
                    ++log_event_count;
                    if (log_event_count % 10 == 0) {
                        const int total_cuts = base.nrows - R.wm.nrows;
                        fprintf(stderr, "           model: %d rows x %d cols, %d cuts\n",
                                base.nrows, base.ncols, total_cuts);
                    }
                }
                next_log_sec = bnb_elapsed + config.log_interval_seconds;
            }
        }

        // Get next node
        if (frontier.empty()) {
            out.frontier_exhausted = true;
            break;
        }
        const int node_id = frontier.pop();
        BranchNodeState branch_node = nodes[static_cast<size_t>(node_id)];

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
        if (branch_node.parent_dual_bound >= R.cutoff() - tol) {
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

        const auto &effective_decisions = augmented_decisions.empty()
                                              ? branch_node.decisions
                                              : augmented_decisions;

        // Decomposition: if all linking columns are decided, solve sub-problems.
        if (decomp_state.enabled && all_linkers_fixed(decomp_state, effective_decisions)) {
            const double remaining = (config.time_limit_seconds > 0)
                ? std::max(1.0, config.time_limit_seconds - R.elapsed())
                : 60.0;

            auto dr = solve_decomposed(base, effective_decisions, decomp_state,
                                       R.cutoff(), config, remaining, verbosity);
            if (dr.solved) {
                // Sub-problems fully explored this subtree
                if (dr.improved && R.try_adopt(dr.combined_obj, dr.combined_solution, "decomposition"))
                    on_new_incumbent(nullptr, 0.0, 0.0);
                continue; // skip LP solve/branch: subtree handled
            }
            // Decomposition failed (infeasible block or no split):
            // fall through to normal LP solve + branching
        }

        // Solve LP with per-node basis warm-start
        lp.apply_decisions(effective_decisions);
        lp.add_cuts(branch_node.cuts);
        if (has_node_basis) {
            lp.set_basis(basis_it->second);
            node_bases.erase(basis_it);
        } else {
            lp.restore_basis();
        }
        LpSolution sol = lp.solve();
        ++R.total_lp_solves;
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
        if (R.obj_is_integral && sol.optimal && std::isfinite(node_dual_bound))
            node_dual_bound = tighten_dual_bound(node_dual_bound, integ_tol);

        const bool dual_improved = sol.optimal &&
            (node_dual_bound > branch_node.parent_dual_bound + tol);

        // Heuristics
        if (processed_nodes == 1 ||
            (config.heuristic_frequency > 0 && processed_nodes % config.heuristic_frequency == 0) ||
            dual_improved) {
            if (R.run_heuristics(sol, branch_node, "", true) &&
                on_new_incumbent(&branch_node, node_dual_bound, node_dual_bound_raw))
                continue;
        }

        // Node-level Lagrangian heuristic
        if (config.lagrangian_frequency > 0.0 && sol.optimal &&
            static_cast<int>(sol.row_dual.size()) >= base.nrows) {
            lagrangian_accumulator += config.lagrangian_frequency;
            if (lagrangian_accumulator >= 1.0) {
                lagrangian_accumulator -= 1.0;

                auto lagr = lagrangian_relaxation_at_node(
                    base, R.cutoff(), sol.row_dual, effective_decisions,
                    50,   // max_iterations (reduced from 500)
                    0.5); // time_limit (reduced from 5.0)

                if (!lagr.best_solution.empty() &&
                    R.try_adopt(lagr.best_heuristic_obj, lagr.best_solution, "node_lagrangian") &&
                    on_new_incumbent(&branch_node, node_dual_bound, node_dual_bound_raw))
                    continue;
            }
        }

        // Diving heuristic (skip in proving phase: frontier shrinking)
        if (config.diving_frequency > 0.0 && sol.optimal && !proving_phase &&
            node_dual_bound < R.cutoff() - tol) {
            diving_accumulator += config.diving_frequency;
            if (diving_accumulator >= 1.0) {
                diving_accumulator -= 1.0;

                auto dive = run_diving_heuristics(
                    lp, base, branch_node, sol, R.cutoff(), integ_tol,
                    config.diving_max_lp_solves);
                R.total_lp_solves += dive.lp_solves;

                if (dive.found &&
                    R.try_adopt(dive.objective, dive.solution, dive.strategy_name) &&
                    on_new_incumbent(&branch_node, node_dual_bound, node_dual_bound_raw))
                    continue;
            }
        }

        // RINS heuristic (skip in proving phase: frontier shrinking)
        if (config.rins_frequency > 0.0 && sol.optimal && !proving_phase &&
            !R.inc.solution.empty() && node_dual_bound < R.cutoff() - tol) {
            rins_accumulator += config.rins_frequency;
            if (rins_accumulator >= 1.0) {
                rins_accumulator -= 1.0;

                // Incumbent restricted to the active columns
                incumbent_active.assign(static_cast<size_t>(base.ncols), 0.0);
                for (int j = 0; j < base.ncols; ++j) {
                    const int orig = base.active_to_original[static_cast<size_t>(j)];
                    if (orig >= 0 && orig < R.ncols_input &&
                        R.inc.solution[static_cast<size_t>(orig)] > 0.5)
                        incumbent_active[static_cast<size_t>(j)] = 1.0;
                }

                const double remaining = (config.time_limit_seconds > 0)
                    ? std::max(1.0, config.time_limit_seconds - R.elapsed())
                    : 30.0;

                auto rr = run_rins(base, sol, incumbent_active, R.cutoff(),
                                   config, std::min(remaining, 30.0), verbosity);
                R.total_lp_solves += rr.sub_lp_solves;

                if (rr.found &&
                    R.try_adopt(rr.objective, rr.solution, "rins") &&
                    on_new_incumbent(&branch_node, node_dual_bound, node_dual_bound_raw))
                    continue;
            }
        }

        // Pruning by bound
        if (node_dual_bound >= R.cutoff() - tol) continue;

        // Integrality check
        if (sol.optimal && is_binary_integral(sol.col_value, base.ncols, integ_tol)) {
            if (R.try_adopt(sol.primal_obj, sol.col_value, "exact_node"))
                on_new_incumbent(nullptr, 0.0, 0.0);
            continue;
        }

        // Compute reduced costs once per node (used by RC fixing and Balas)
        std::vector<double> rcosts;
        const bool need_rcosts = sol.optimal &&
            ((std::isfinite(R.cutoff()) && node_dual_bound < R.cutoff() - tol) ||
             (config.balas_enabled && force_aggressive_branching));
        if (need_rcosts)
            rcosts = compute_reduced_costs(base.obj, sol.row_dual, base, base.ncols);

        // Node-level reduced cost fixing: fix variables that cannot improve
        std::vector<BranchDecision> rc_fixings;
        if (!rcosts.empty() && std::isfinite(R.cutoff()) && node_dual_bound < R.cutoff() - tol) {
            const double node_gap = R.cutoff() - node_dual_bound;
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
                                            R.cutoff(), base, base.ncols, config.balas_max_branches,
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
                if (verbosity >= 3)
                    fprintf(stderr, "  [%8.3fs] Aggressive BR1: %d children from %d sets\n",
                            R.elapsed(), enqueued, static_cast<int>(br.sets.size()));
            }
            force_aggressive_branching = false;
        }

        if (!used_balas) {
            int branch_var = -1;

            // Priority: branch on unfixed linking columns first
            if (decomp_state.enabled)
                branch_var = pick_linking_branch_var(decomp_state, effective_decisions, fractional);

            if (branch_var < 0) {
                if (use_reliability) {
                    int sb_lps = 0;
                    if (!adj_valid) {
                        cached_adj = build_adjacency(base);
                        adj_valid = true;
                    }
                    branch_var = reliability_branch_select(
                        pc_state, lp, base, branch_node, sol, fractional,
                        R.cutoff(), config.reliability_eta, config.reliability_max_sb,
                        integ_tol, verbosity, sb_lps, &cached_adj);
                    R.total_lp_solves += sb_lps;
                } else {
                    branch_var = selector->select(sol.col_value, base.obj, fractional);
                }
            }
            if (branch_var < 0) continue;

            for (int value = 0; value <= 1; ++value) {
                BranchNodeState child;
                if (!append_decision_if_consistent(branch_node, branch_var, value, &child)) continue;
                for (const auto &rc : rc_fixings)
                    child.decisions.push_back(rc);
                if (is_node_provably_infeasible(child, base)) continue;
                child.parent_dual_bound = node_dual_bound;
                child.parent_dual_bound_raw = node_dual_bound_raw;
                const int child_id = static_cast<int>(nodes.size());
                nodes.push_back(std::move(child));
                node_bases[child_id] = node_basis;
                frontier.push(child_id);
            }
        }

        // Update global dual bound from heap top: O(1) instead of O(|frontier|)
        if (!frontier.empty()) {
            R.global_dual_bound = frontier.top_bound();
            R.global_dual_bound_raw = frontier.min_raw_bound;
        }

        // Stagnation control
        if (gap_stagnation_window > 0 && std::isfinite(R.inc.obj)) {
            const double current_gap = compute_mip_gap(R.inc.obj, R.global_dual_bound + R.fixed_cost);
            if (std::isfinite(current_gap) && current_gap < best_mip_gap_seen - 1e-8) {
                best_mip_gap_seen = current_gap;
                node_at_last_gap_improvement = processed_nodes;
            }

            if (processed_nodes - node_at_last_gap_improvement >= gap_stagnation_window) {
                node_at_last_gap_improvement = processed_nodes;

                // Dynamic decomposition activation after ~150s of BnB
                if (!decomp_state.enabled && decomp_dynamic_eligible) {
                    if (R.elapsed() > 150.0) {
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
                    if (verbosity >= 2)
                        fprintf(stderr, "  [%8.3fs] Proving phase: frontier shrinking, "
                                        "disabling expensive heuristics\n", R.elapsed());
                }
                if (proving_phase && frontier_shrink_streak <= -3) {
                    proving_phase = false;
                    if (verbosity >= 2)
                        fprintf(stderr, "  [%8.3fs] Exiting proving phase: frontier growing\n",
                                R.elapsed());
                }
                const bool frontier_shrinking = (frontier_shrink_streak >= 2);
                frontier_at_last_stagnation = frontier.size();

                // Mid-BnB cuts (skip if frontier is naturally shrinking)
                cut_accumulator += config.mid_bnb_cut_frequency;
                if (cut_accumulator >= 1.0 && !frontier_shrinking &&
                    config.cuts_enabled && sol.optimal) {
                    cut_accumulator -= 1.0;

                    // Save pre-cut state for rollback
                    const double pre_cut_dual = R.global_dual_bound;
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
                                                          base, base.ncols, integ_tol, R.cutoff());
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
                        ++R.total_lp_solves;
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
                        if (R.obj_is_integral && std::isfinite(post_cut_dual))
                            post_cut_dual = tighten_dual_bound(post_cut_dual, integ_tol);
                        if (post_cut_dual <= pre_cut_dual + tol) {
                            // Cuts didn't help: roll back
                            base.nrows = pre_cut_nrows;
                            base.nnz = pre_cut_nnz;
                            base.csr_inds = pre_cut_csr_inds;
                            base.csr_offs = pre_cut_csr_offs;
                            base.csr_vals = pre_cut_csr_vals;
                            base.rhs = pre_cut_rhs;
                            base.base_cuts = pre_cut_base_cuts;
                            lp.rebuild_model(base);
                            if (verbosity >= 3)
                                fprintf(stderr, "  [%8.3fs] Stagnation: %d mid-BnB cuts undone (no dual improvement)\n",
                                        R.elapsed(), total_cuts);
                        } else if (verbosity >= 3) {
                            fprintf(stderr, "  [%8.3fs] Stagnation: added %d mid-BnB cuts (dual %.8f -> %.8f)\n",
                                    R.elapsed(), total_cuts, pre_cut_dual, post_cut_dual);
                        }
                    }
                }

                // Aggressive Balas branching (independent of cuts)
                balas_accumulator += config.aggressive_balas_frequency;
                if (balas_accumulator >= 1.0 && config.balas_enabled) {
                    balas_accumulator -= 1.0;
                    force_aggressive_branching = true;
                    if (verbosity >= 3)
                        fprintf(stderr, "  [%8.3fs] Stagnation: aggressive Balas branching\n", R.elapsed());
                }
            }
        }
    }

    // Final dual bound from the remaining frontier (reduced space)
    {
        double new_bound = kInf;
        double new_bound_raw = kInf;
        for (const int idx : frontier.heap) {
            new_bound = std::min(new_bound, nodes[static_cast<size_t>(idx)].parent_dual_bound);
            new_bound_raw = std::min(new_bound_raw, nodes[static_cast<size_t>(idx)].parent_dual_bound_raw);
        }
        if (std::isfinite(new_bound)) R.global_dual_bound = new_bound;
        else if (frontier.empty() && std::isfinite(R.inc.obj)) R.global_dual_bound = R.cutoff();
        if (std::isfinite(new_bound_raw)) R.global_dual_bound_raw = new_bound_raw;
        else if (frontier.empty() && std::isfinite(R.inc.obj)) R.global_dual_bound_raw = R.cutoff();
    }

    out.processed_nodes = processed_nodes;
    return out;
}

SolverResult make_result(const SolverRun &R, const BnbOutcome &out) {
    SolverResult result;
    result.wall_time = R.elapsed();
    result.nodes_processed = out.processed_nodes;
    result.lp_solves = R.total_lp_solves;

    // Dual bound in the original space
    const double dual_total = R.global_dual_bound + R.fixed_cost;

    if (std::isfinite(R.inc.obj)) {
        result.primal_obj = R.inc.obj;
        result.solution = R.inc.solution;

        if ((out.frontier_exhausted || out.gap_tolerance_reached) &&
            !out.hard_time_limit_reached &&
            out.processed_nodes < R.config.max_nodes) {
            result.dual_obj = R.inc.obj;
            result.mip_gap = 0.0;
            result.status = "Optimal";
        } else {
            result.dual_obj = dual_total;
            result.mip_gap = compute_mip_gap(R.inc.obj, dual_total);
            result.status = out.hard_time_limit_reached ? "TimeLimit" : "NodeLimit";
        }
    } else {
        result.primal_obj = kInf;
        result.dual_obj = dual_total;
        result.mip_gap = kInf;
        result.status = "NoSolution";
    }
    return result;
}

} // namespace

SolverResult solve(const ScpInstance &instance, const SolverConfig &config) {
    SolverRun R(instance, config);
    if (R.obj_is_integral && R.verbosity >= 2)
        fprintf(stderr, "Objective coefficients are integral; enabling dual bound tightening\n");

    run_greedy(R);
    run_root_reductions(R);
    if (R.wm.nrows == 0)
        return trivial_result(R);

    solve_root_lp(R);
    run_post_lp_phase(R);
    prepare_bnb_model(R);

    const BnbOutcome out = run_branch_and_bound(R);
    return make_result(R, out);
}

} // namespace scpsol
