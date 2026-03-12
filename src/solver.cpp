#include "solver.h"
#include "balas.h"
#include "bnb.h"
#include "cuts.h"
#include "heuristics.h"
#include "lp.h"
#include "preprocessor.h"
#include "reliability.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <deque>
#include <limits>
#include <memory>
#include <string>

namespace scpsol {

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

static void prune_frontier(
    std::deque<int> &frontier,
    const std::vector<BranchNodeState> &nodes,
    double best_obj,
    double tol,
    int verbosity) {
    const double prune_bound = best_obj - tol;
    size_t before = frontier.size();
    std::deque<int> surviving;
    for (const int idx : frontier) {
        if (nodes[static_cast<size_t>(idx)].parent_dual_bound < prune_bound)
            surviving.push_back(idx);
    }
    frontier.swap(surviving);
    if (frontier.size() < before && verbosity >= 3) {
        fprintf(stderr, "  Frontier pruned: %zu -> %zu nodes\n", before, frontier.size());
    }
}

static int mid_bnb_column_removal(
    BaseRelaxationModel &base,
    double best_obj,
    double tol,
    std::deque<int> &frontier,
    std::vector<BranchNodeState> &nodes,
    int verbosity,
    PseudocostState *pc_state = nullptr) {
    ModelReductionResult reduction = reduce_base_model(base, best_obj, tol);
    if (reduction.columns_removed <= 0) return 0;
    if (verbosity >= 3)
        fprintf(stderr, "  Mid-BnB reduction: %d cols removed, %d remaining\n",
                reduction.columns_removed, base.ncols);
    std::deque<int> surviving;
    for (const int idx : frontier) {
        if (remap_branch_node(nodes[static_cast<size_t>(idx)], reduction.old_to_new))
            surviving.push_back(idx);
    }
    frontier.swap(surviving);
    if (pc_state) pc_state->remap(reduction.old_to_new, base.ncols);
    return reduction.columns_removed;
}

static int mid_bnb_budget_pruning(
    BaseRelaxationModel &base,
    double best_obj,
    double tol,
    double preprocess_time_limit,
    std::deque<int> &frontier,
    std::vector<BranchNodeState> &nodes,
    int verbosity,
    PseudocostState *pc_state = nullptr) {
    ModelReductionResult reduction = reduce_base_model_budget_pruning(
        base, best_obj, tol, preprocess_time_limit);
    if (reduction.columns_removed <= 0) return 0;
    if (verbosity >= 3)
        fprintf(stderr, "  Mid-BnB budget pruning: %d cols removed, %d remaining\n",
                reduction.columns_removed, base.ncols);
    std::deque<int> surviving;
    for (const int idx : frontier) {
        if (remap_branch_node(nodes[static_cast<size_t>(idx)], reduction.old_to_new))
            surviving.push_back(idx);
    }
    frontier.swap(surviving);
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
    // Phase 2: Cost + budget + dominance reduction
    // ================================================================
    {
        const int cols_before = ncols;
        apply_cost_reduction(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                             best_obj, tol, verbosity);
        apply_budget_pruning_preprocess(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                                        best_obj, tol, config.preprocess_time_limit, verbosity);
        apply_dominance_reduction(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                                  "cost_driven", tol, config.preprocess_time_limit, verbosity);
        apply_dominance_reduction(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                                  config.preprocess_rules, tol, config.preprocess_time_limit, verbosity);
        // Row reduction: essential columns, row domination, probing
        auto rr = row_reduce(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                             tol, config.preprocess_time_limit, verbosity);
        if (rr.cols_fixed > 0 || rr.rows_removed > 0) {
            fixed_preprocess_cost += rr.fixed_cost;
            for (int c : rr.fixed_original_cols) fixed_original_cols.push_back(c);
            best_obj -= rr.fixed_cost;
            if (verbosity >= 2)
                fprintf(stderr, "  Row reduction: %d rows removed, %d cols fixed (cost %.6g)\n",
                        rr.rows_removed, rr.cols_fixed, rr.fixed_cost);
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
    };
    build_base();

    LpSolver lp;
    lp.build_model(base);
    int total_lp_solves = 0;

    LpSolution root_sol = lp.solve();
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
    // Phase 4-5: Post-LP reduction
    // ================================================================
    {
        const int cols_before = ncols;
        // Update working arrays from base
        csr_inds = base.csr_inds; csr_offs = base.csr_offs; csr_vals = base.csr_vals;
        obj = base.obj; active_to_input = base.active_to_original;
        ncols = base.ncols;

        // Reduced-cost fixing: remove columns with reduced cost > gap
        if (root_sol.solved && root_sol.optimal && std::isfinite(best_obj)) {
            const double gap = best_obj - root_sol.dual_obj;
            if (gap > tol) {
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
                if (rc_removed > 0) {
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
                        fprintf(stderr, "  Reduced-cost fixing: %d cols removed, %d remaining\n",
                                rc_removed, nc);
                    ncols = nc;
                    csr_inds = std::move(new_inds_rc);
                    csr_offs = std::move(new_offs_rc);
                    csr_vals = std::move(new_vals_rc);
                    obj = std::move(new_obj_rc);
                    active_to_input = std::move(new_active_input);
                }
            }
        }

        apply_cost_reduction(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                             best_obj, tol, verbosity);
        apply_budget_pruning_preprocess(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                                        best_obj, tol, config.preprocess_time_limit, verbosity);
        apply_dominance_reduction(nrows, ncols, csr_inds, csr_offs, csr_vals, obj, active_to_input,
                                  config.preprocess_rules, tol, config.preprocess_time_limit, verbosity);
        // Second round of row reduction after LP-based column removal
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
                if (verbosity >= 2)
                    fprintf(stderr, "  Post-LP row reduction: %d rows removed, %d cols fixed (cost %.6g)\n",
                            rr2.rows_removed, rr2.cols_fixed, rr2.fixed_cost);
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

    if (verbosity >= 2)
        fprintf(stderr, "BnB base model: %d rows x %d cols, nnz=%d\n",
                base.nrows, base.ncols, base.nnz);

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

    std::vector<BranchNodeState> nodes;
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

    std::deque<int> frontier;
    frontier.push_back(0);

    int processed_nodes = 0;
    bool gap_tolerance_reached = false;
    bool hard_time_limit_reached = false;
    bool frontier_exhausted = false;

    const int gap_stagnation_window = config.gap_stagnation_window;
    double best_mip_gap_seen = std::numeric_limits<double>::infinity();
    int node_at_last_gap_improvement = 0;

    double cut_accumulator = 0.0;
    double balas_accumulator = 0.0;
    size_t frontier_at_last_stagnation = 0;
    int frontier_shrink_streak = 0;
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
                fprintf(stderr, "BnB hard time limit reached (%.1fs)\n", elapsed);
            break;
        }

        // Gap check
        if (std::isfinite(best_obj) && std::isfinite(global_dual_bound)) {
            const double gap = compute_mip_gap(best_obj, global_dual_bound);
            if (std::isfinite(gap) && gap <= config.mip_gap_tol) {
                gap_tolerance_reached = true;
                if (verbosity >= 1)
                    fprintf(stderr, "MIP gap %.6f%% within tolerance; optimal\n", gap * 100.0);
                break;
            }
        }

        // Log
        if (config.log_interval_seconds > 0.0) {
            const double bnb_elapsed = std::chrono::duration<double>(now - bnb_start).count();
            if (bnb_elapsed >= next_log_sec) {
                if (verbosity >= 2) {
                    const double gap = compute_mip_gap(best_obj, global_dual_bound_raw);
                    fprintf(stderr, "  nodes=%4d frontier=%4zu lp=%5d incumbent=%.12g dual=%.12g gap=",
                            processed_nodes, frontier.size(), total_lp_solves,
                            best_obj, global_dual_bound_raw);
                    if (std::isfinite(gap)) fprintf(stderr, "%.4f%%\n", gap * 100.0);
                    else fprintf(stderr, "inf\n");
                    ++log_event_count;
                    if (log_event_count % 10 == 0) {
                        int total_cuts = base.nrows - nrows;
                        fprintf(stderr, "    model: %d rows x %d cols, %d cuts\n",
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
        const int node_id = frontier.front();
        frontier.pop_front();
        const BranchNodeState branch_node = nodes[static_cast<size_t>(node_id)];

        // Bound check
        if (branch_node.parent_dual_bound >= best_obj - tol) continue;

        // Solve LP
        lp.apply_decisions(branch_node.decisions);
        lp.add_cuts(branch_node.cuts);
        lp.restore_basis();
        LpSolution sol = lp.solve();
        ++total_lp_solves;
        lp.save_basis();
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
                    if (verbosity >= 2)
                        fprintf(stderr, "  New incumbent: %.12g (from %s)\n", best_obj, hr.name.c_str());
                    break;
                }
            }
            if (incumbent_improved) {
                node_at_last_gap_improvement = processed_nodes;
                prune_frontier(frontier, nodes, best_obj, tol, verbosity);
                if (mid_bnb_column_removal(base, best_obj, tol, frontier, nodes, verbosity, &pc_state) > 0)
                    lp.rebuild_model_keep_basis(base);
                if (mid_bnb_budget_pruning(base, best_obj, tol, config.preprocess_time_limit, frontier, nodes, verbosity, &pc_state) > 0)
                    lp.rebuild_model_keep_basis(base);
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
                if (verbosity >= 2)
                    fprintf(stderr, "  New incumbent: %.12g (exact)\n", best_obj);
                prune_frontier(frontier, nodes, best_obj, tol, verbosity);
                if (mid_bnb_column_removal(base, best_obj, tol, frontier, nodes, verbosity, &pc_state) > 0)
                    lp.rebuild_model_keep_basis(base);
                if (mid_bnb_budget_pruning(base, best_obj, tol, config.preprocess_time_limit, frontier, nodes, verbosity, &pc_state) > 0)
                    lp.rebuild_model_keep_basis(base);
            }
            continue;
        }

        // Branch
        auto fractional = collect_fractional_candidates(sol.col_value, base.ncols, integ_tol);
        if (fractional.empty()) continue;

        bool used_balas = false;
        if (config.balas_enabled && force_aggressive_branching) {
            auto rcosts = compute_reduced_costs(base.obj, sol.row_dual, base, base.ncols);
            auto br = balas_branch_generate(sol.col_value, sol.row_dual, rcosts,
                                            best_obj, base, base.ncols, config.balas_max_branches,
                                            integ_tol, true);
            if (br.use_br1) {
                auto children = balas_br1_children(branch_node, br.sets, node_dual_bound, node_dual_bound_raw);
                int enqueued = 0;
                for (auto &child : children) {
                    if (!is_node_provably_infeasible(child, base)) {
                        nodes.push_back(std::move(child));
                        frontier.push_back(static_cast<int>(nodes.size()) - 1);
                        ++enqueued;
                    }
                }
                used_balas = (enqueued > 0);
                if (verbosity >= 3) {
                    fprintf(stderr, "  Aggressive BR1: %d children from %d sets\n",
                            enqueued, static_cast<int>(br.sets.size()));
                }
            }
            force_aggressive_branching = false;
        }

        if (!used_balas) {
            int branch_var;
            if (use_reliability) {
                int sb_lps = 0;
                branch_var = reliability_branch_select(
                    pc_state, lp, base, branch_node, sol, fractional,
                    best_obj, config.reliability_eta, config.reliability_max_sb,
                    integ_tol, verbosity, sb_lps);
                total_lp_solves += sb_lps;
            } else {
                branch_var = selector->select(sol.col_value, base.obj, fractional);
            }
            if (branch_var < 0) continue;

            BranchNodeState child_zero;
            if (append_decision_if_consistent(branch_node, branch_var, 0, &child_zero) &&
                !is_node_provably_infeasible(child_zero, base)) {
                child_zero.parent_dual_bound = node_dual_bound;
                child_zero.parent_dual_bound_raw = node_dual_bound_raw;
                nodes.push_back(child_zero);
                frontier.push_back(static_cast<int>(nodes.size()) - 1);
            }

            BranchNodeState child_one;
            if (append_decision_if_consistent(branch_node, branch_var, 1, &child_one) &&
                !is_node_provably_infeasible(child_one, base)) {
                child_one.parent_dual_bound = node_dual_bound;
                child_one.parent_dual_bound_raw = node_dual_bound_raw;
                nodes.push_back(child_one);
                frontier.push_back(static_cast<int>(nodes.size()) - 1);
            }
        }

        // Update global dual bound from frontier
        {
            double new_bound = std::numeric_limits<double>::infinity();
            double new_bound_raw = std::numeric_limits<double>::infinity();
            for (const int idx : frontier) {
                new_bound = std::min(new_bound, nodes[static_cast<size_t>(idx)].parent_dual_bound);
                new_bound_raw = std::min(new_bound_raw, nodes[static_cast<size_t>(idx)].parent_dual_bound_raw);
            }
            if (std::isfinite(new_bound)) global_dual_bound = new_bound;
            if (std::isfinite(new_bound_raw)) global_dual_bound_raw = new_bound_raw;
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

                // Track frontier trend
                bool frontier_shrinking = false;
                if (frontier_at_last_stagnation > 0) {
                    if (frontier.size() < frontier_at_last_stagnation)
                        ++frontier_shrink_streak;
                    else
                        frontier_shrink_streak = 0;
                    frontier_shrinking = (frontier_shrink_streak >= 2);
                }
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
                            if (verbosity >= 3)
                                fprintf(stderr, "  Stagnation: %d mid-BnB cuts undone (no dual improvement)\n",
                                        total_cuts);
                        } else {
                            if (verbosity >= 3)
                                fprintf(stderr, "  Stagnation: added %d mid-BnB cuts (dual %.6g -> %.6g)\n",
                                        total_cuts, pre_cut_dual, post_cut_dual);
                        }
                    }
                }

                // Aggressive Balas branching (independent of cuts)
                balas_accumulator += config.aggressive_balas_frequency;
                if (balas_accumulator >= 1.0 && config.balas_enabled) {
                    balas_accumulator -= 1.0;
                    force_aggressive_branching = true;
                    if (verbosity >= 3)
                        fprintf(stderr, "  Stagnation: aggressive Balas branching\n");
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
        for (const int idx : frontier) {
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
