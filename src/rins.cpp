#include "rins.h"
#include "solver.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

namespace scpsol {

RinsResult run_rins(
    const BaseRelaxationModel &base,
    const LpSolution &lp_sol,
    const std::vector<double> &incumbent_active,
    double incumbent_obj,
    const SolverConfig &config,
    double time_budget,
    int verbosity) {

    RinsResult result;
    const int ncols = base.ncols;
    const int nrows = base.nrows;

    if (ncols <= 0 || nrows <= 0) return result;
    if (static_cast<int>(lp_sol.col_value.size()) < ncols) return result;
    if (static_cast<int>(incumbent_active.size()) < ncols) return result;

    const double tol = 0.1; // agreement threshold

    // Classify columns: agree-on-1, agree-on-0, or free (disagree)
    std::vector<int8_t> status(static_cast<size_t>(ncols), 0);
    // +1 = fixed to 1 (both agree), -1 = fixed to 0 (both agree), 0 = free
    int fixed_one = 0, fixed_zero = 0, free_count = 0;
    double fixed_cost = 0.0;

    for (int j = 0; j < ncols; ++j) {
        const double lp_val = lp_sol.col_value[static_cast<size_t>(j)];
        const double inc_val = incumbent_active[static_cast<size_t>(j)];

        if (lp_val > 1.0 - tol && inc_val > 0.5) {
            // Both ~1: fix to 1
            status[static_cast<size_t>(j)] = 1;
            fixed_cost += base.obj[static_cast<size_t>(j)];
            ++fixed_one;
        } else if (lp_val < tol && inc_val < 0.5) {
            // Both ~0: fix to 0
            status[static_cast<size_t>(j)] = -1;
            ++fixed_zero;
        } else {
            // Disagree: free variable for sub-problem
            ++free_count;
        }
    }

    // Skip if too few variables are fixed (RINS is only useful when
    // most variables agree) or if all are fixed, or if the sub-problem
    // is too large relative to the original
    if (free_count == 0 || fixed_one + fixed_zero < ncols / 3 ||
        free_count > ncols / 2)
        return result;

    // Determine which rows are already covered by fixed-to-1 columns
    std::vector<bool> row_covered(static_cast<size_t>(nrows), false);
    for (int i = 0; i < nrows; ++i) {
        for (int k = base.csr_offs[static_cast<size_t>(i)];
             k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int c = base.csr_inds[static_cast<size_t>(k)];
            if (c >= 0 && c < ncols && status[static_cast<size_t>(c)] == 1) {
                row_covered[static_cast<size_t>(i)] = true;
                break;
            }
        }
    }

    // Collect uncovered rows and free columns that cover them
    std::vector<int> sub_rows, sub_cols;
    for (int i = 0; i < nrows; ++i)
        if (!row_covered[static_cast<size_t>(i)])
            sub_rows.push_back(i);

    if (sub_rows.empty()) {
        // All rows covered by fixed columns — the fixed solution is feasible
        result.found = true;
        result.objective = fixed_cost;
        result.solution.assign(static_cast<size_t>(ncols), 0.0);
        for (int j = 0; j < ncols; ++j)
            if (status[static_cast<size_t>(j)] == 1)
                result.solution[static_cast<size_t>(j)] = 1.0;
        return result;
    }

    // Identify which free columns can cover uncovered rows
    std::vector<bool> col_useful(static_cast<size_t>(ncols), false);
    for (int i : sub_rows) {
        for (int k = base.csr_offs[static_cast<size_t>(i)];
             k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int c = base.csr_inds[static_cast<size_t>(k)];
            if (c >= 0 && c < ncols && status[static_cast<size_t>(c)] == 0)
                col_useful[static_cast<size_t>(c)] = true;
        }
    }
    for (int j = 0; j < ncols; ++j)
        if (col_useful[static_cast<size_t>(j)])
            sub_cols.push_back(j);

    if (sub_cols.empty())
        return result; // infeasible sub-problem

    // Build ScpInstance for the sub-problem
    ScpInstance sub;
    sub.nrows = static_cast<int>(sub_rows.size());
    sub.ncols = static_cast<int>(sub_cols.size());

    // Column mapping
    std::vector<int> sub_to_parent(sub_cols.begin(), sub_cols.end());
    std::vector<int> parent_to_sub(static_cast<size_t>(ncols), -1);
    for (int jj = 0; jj < sub.ncols; ++jj)
        parent_to_sub[static_cast<size_t>(sub_to_parent[static_cast<size_t>(jj)])] = jj;

    // Costs
    sub.costs.resize(static_cast<size_t>(sub.ncols));
    for (int jj = 0; jj < sub.ncols; ++jj)
        sub.costs[static_cast<size_t>(jj)] =
            base.obj[static_cast<size_t>(sub_to_parent[static_cast<size_t>(jj)])];

    // Build CSR
    sub.csr_offsets.resize(static_cast<size_t>(sub.nrows) + 1, 0);
    for (int ii = 0; ii < sub.nrows; ++ii) {
        const int parent_row = sub_rows[static_cast<size_t>(ii)];
        for (int k = base.csr_offs[static_cast<size_t>(parent_row)];
             k < base.csr_offs[static_cast<size_t>(parent_row) + 1]; ++k) {
            const int pc = base.csr_inds[static_cast<size_t>(k)];
            if (pc >= 0 && pc < ncols &&
                parent_to_sub[static_cast<size_t>(pc)] >= 0) {
                sub.csr_indices.push_back(parent_to_sub[static_cast<size_t>(pc)]);
                sub.csr_values.push_back(base.csr_vals[static_cast<size_t>(k)]);
            }
        }
        sub.csr_offsets[static_cast<size_t>(ii) + 1] =
            static_cast<int>(sub.csr_indices.size());
    }

    if (verbosity >= 3)
        fprintf(stderr, "    RINS: %d fixed-1, %d fixed-0, %d free -> sub %dx%d\n",
                fixed_one, fixed_zero, free_count, sub.nrows, sub.ncols);

    // Budget for sub-problem: remaining incumbent minus fixed cost
    double sub_budget = incumbent_obj - fixed_cost;
    if (sub_budget <= 0)
        return result; // fixed cost already exceeds incumbent

    // Solve sub-problem
    SolverConfig sub_config = config;
    sub_config.verbosity = std::max(0, verbosity - 2);
    sub_config.time_limit_seconds = std::min(time_budget, 5.0);
    sub_config.max_nodes = std::min(config.max_nodes, 5000);
    sub_config.decomposition_mode = "off";
    sub_config.diving_frequency = 0.0;
    sub_config.lagrangian_frequency = 0.0;
    // Keep RINS disabled in sub-problem to avoid recursion
    // (controlled by caller not passing rins_frequency to sub-solver)

    SolverResult sub_result = solve(sub, sub_config);
    result.sub_nodes = sub_result.nodes_processed;
    result.sub_lp_solves = sub_result.lp_solves;

    if (sub_result.status == "Infeasible" || !std::isfinite(sub_result.primal_obj))
        return result;

    double combined_obj = fixed_cost + sub_result.primal_obj;
    if (combined_obj >= incumbent_obj - 1e-6)
        return result; // no improvement

    // Build combined solution in parent model space
    result.found = true;
    result.objective = combined_obj;
    result.solution.assign(static_cast<size_t>(ncols), 0.0);

    // Fixed-to-1 columns
    for (int j = 0; j < ncols; ++j)
        if (status[static_cast<size_t>(j)] == 1)
            result.solution[static_cast<size_t>(j)] = 1.0;

    // Sub-problem solution
    for (int jj = 0; jj < sub.ncols; ++jj) {
        if (static_cast<size_t>(jj) < sub_result.solution.size() &&
            sub_result.solution[static_cast<size_t>(jj)] > 0.5) {
            result.solution[static_cast<size_t>(
                sub_to_parent[static_cast<size_t>(jj)])] = 1.0;
        }
    }

    if (verbosity >= 3)
        fprintf(stderr, "    RINS: sub solved (%s) obj=%.0f, combined=%.0f (incumbent=%.0f)\n",
                sub_result.status.c_str(), sub_result.primal_obj,
                combined_obj, incumbent_obj);

    return result;
}

} // namespace scpsol
