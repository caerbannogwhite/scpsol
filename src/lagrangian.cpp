#include "lagrangian.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>

namespace scpsol {

LagrangianResult lagrangian_relaxation(
    const BaseRelaxationModel &base,
    double incumbent_obj,
    const std::vector<double> &init_multipliers,
    int max_iterations,
    double time_limit,
    int verbosity) {

    using Clock = std::chrono::steady_clock;
    const auto start = Clock::now();

    const int nrows = base.nrows;
    const int ncols = base.ncols;
    LagrangianResult result;

    if (nrows <= 0 || ncols <= 0) return result;

    // Initialize multipliers from LP duals (clamped >= 0)
    std::vector<double> u(static_cast<size_t>(nrows), 0.0);
    if (static_cast<int>(init_multipliers.size()) >= nrows) {
        for (int i = 0; i < nrows; ++i)
            u[static_cast<size_t>(i)] = std::max(0.0, init_multipliers[static_cast<size_t>(i)]);
    }

    const auto &cols_to_rows = base.cols_to_rows;
    if (static_cast<int>(cols_to_rows.size()) < ncols) return result;

    double lambda = 2.0;
    int no_improve_count = 0;
    const int halve_interval = 30;
    double best_bound = -std::numeric_limits<double>::infinity();
    double ub = incumbent_obj; // best known upper bound

    std::vector<double> x(static_cast<size_t>(ncols));         // Lagrangian solution
    std::vector<double> rc(static_cast<size_t>(ncols));        // reduced costs
    std::vector<double> subgrad(static_cast<size_t>(nrows));   // subgradient
    std::vector<double> coverage(static_cast<size_t>(nrows));  // row coverage for heuristic

    for (int iter = 0; iter < max_iterations; ++iter) {
        const double elapsed = std::chrono::duration<double>(Clock::now() - start).count();
        if (elapsed > time_limit) break;

        // --- Lagrangian subproblem ---
        // L(u) = min (c - A'u)'x + u'rhs  s.t. 0 <= x <= 1
        // For each j: rc_j = c_j - sum_i(u_i * a_ij)
        //   x_j = 1 if rc_j < 0, else x_j = 0
        double lagrangian_bound = 0.0;

        // Sum of u_i * rhs_i
        for (int i = 0; i < nrows; ++i)
            lagrangian_bound += u[static_cast<size_t>(i)] * base.rhs[static_cast<size_t>(i)];

        for (int j = 0; j < ncols; ++j) {
            double r = base.obj[static_cast<size_t>(j)];
            for (const auto &entry : cols_to_rows[static_cast<size_t>(j)])
                r -= u[static_cast<size_t>(entry.row)] * entry.val;
            rc[static_cast<size_t>(j)] = r;

            if (r < 0.0) {
                x[static_cast<size_t>(j)] = 1.0;
                lagrangian_bound += r;
            } else {
                x[static_cast<size_t>(j)] = 0.0;
            }
        }

        // Update best bound
        if (lagrangian_bound > best_bound) {
            best_bound = lagrangian_bound;
            no_improve_count = 0;
        } else {
            ++no_improve_count;
        }

        if (no_improve_count >= halve_interval) {
            lambda *= 0.5;
            no_improve_count = 0;
            if (lambda < 1e-6) break;
        }

        // --- Lagrangian heuristic ---
        // Start from Lagrangian solution, repair feasibility, drop redundant columns
        std::fill(coverage.begin(), coverage.end(), 0.0);
        std::vector<double> heur_sol(static_cast<size_t>(ncols), 0.0);

        // Include all columns with x_j = 1 from Lagrangian
        for (int j = 0; j < ncols; ++j) {
            if (x[static_cast<size_t>(j)] > 0.5) {
                heur_sol[static_cast<size_t>(j)] = 1.0;
                for (const auto &entry : cols_to_rows[static_cast<size_t>(j)])
                    coverage[static_cast<size_t>(entry.row)] += entry.val;
            }
        }

        // Greedy repair: cover uncovered rows using Lagrangian reduced costs
        // Prefer columns with smallest (most negative) reduced cost per unit coverage
        while (true) {
            int uncovered = -1;
            for (int i = 0; i < nrows; ++i) {
                if (coverage[static_cast<size_t>(i)] + 1e-8 < base.rhs[static_cast<size_t>(i)]) {
                    uncovered = i;
                    break;
                }
            }
            if (uncovered < 0) break;

            int best_col = -1;
            double best_score = -std::numeric_limits<double>::infinity();
            for (int j = 0; j < ncols; ++j) {
                if (heur_sol[static_cast<size_t>(j)] > 0.5) continue;
                double uncov_gain = 0.0;
                for (const auto &entry : cols_to_rows[static_cast<size_t>(j)]) {
                    if (coverage[static_cast<size_t>(entry.row)] + 1e-8 <
                        base.rhs[static_cast<size_t>(entry.row)] && entry.val > 0.0)
                        uncov_gain += entry.val;
                }
                if (uncov_gain <= 0.0) continue;
                // Score: coverage gain / modified cost (use Lagrangian rc as tiebreaker)
                const double cost = std::max(1e-9, base.obj[static_cast<size_t>(j)]);
                const double score = uncov_gain / cost;
                if (score > best_score) {
                    best_score = score;
                    best_col = j;
                }
            }

            if (best_col < 0) break; // infeasible — shouldn't happen
            heur_sol[static_cast<size_t>(best_col)] = 1.0;
            for (const auto &entry : cols_to_rows[static_cast<size_t>(best_col)])
                coverage[static_cast<size_t>(entry.row)] += entry.val;
        }

        // Drop redundant columns (most expensive first)
        std::vector<int> selected;
        for (int j = 0; j < ncols; ++j) {
            if (heur_sol[static_cast<size_t>(j)] > 0.5)
                selected.push_back(j);
        }
        std::sort(selected.begin(), selected.end(),
                  [&](int a, int b) {
                      return base.obj[static_cast<size_t>(a)] > base.obj[static_cast<size_t>(b)];
                  });

        for (int col : selected) {
            // Try removing this column
            bool can_remove = true;
            for (const auto &entry : cols_to_rows[static_cast<size_t>(col)]) {
                if (coverage[static_cast<size_t>(entry.row)] - entry.val + 1e-8 <
                    base.rhs[static_cast<size_t>(entry.row)]) {
                    can_remove = false;
                    break;
                }
            }
            if (can_remove) {
                heur_sol[static_cast<size_t>(col)] = 0.0;
                for (const auto &entry : cols_to_rows[static_cast<size_t>(col)])
                    coverage[static_cast<size_t>(entry.row)] -= entry.val;
            }
        }

        // Verify feasibility against the CSR matrix directly
        bool feasible = true;
        for (int i = 0; i < nrows; ++i) {
            double row_sum = 0.0;
            for (int k = base.csr_offs[static_cast<size_t>(i)];
                 k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
                const int c = base.csr_inds[static_cast<size_t>(k)];
                if (c >= 0 && c < ncols && heur_sol[static_cast<size_t>(c)] > 0.5)
                    row_sum += base.csr_vals[static_cast<size_t>(k)];
            }
            if (row_sum + 1e-8 < base.rhs[static_cast<size_t>(i)]) {
                feasible = false;
                break;
            }
        }

        if (!feasible) continue;

        // Compute heuristic objective
        double heur_obj = 0.0;
        for (int j = 0; j < ncols; ++j) {
            if (heur_sol[static_cast<size_t>(j)] > 0.5)
                heur_obj += base.obj[static_cast<size_t>(j)];
        }

        if (heur_obj < result.best_heuristic_obj) {
            result.best_heuristic_obj = heur_obj;
            result.best_solution = heur_sol;
            if (heur_obj < ub)
                ub = heur_obj;
        }

        // --- Subgradient update ---
        // s_i = rhs_i - sum_j(a_ij * x_j)
        std::fill(subgrad.begin(), subgrad.end(), 0.0);
        for (int i = 0; i < nrows; ++i)
            subgrad[static_cast<size_t>(i)] = base.rhs[static_cast<size_t>(i)];

        for (int j = 0; j < ncols; ++j) {
            if (x[static_cast<size_t>(j)] < 0.5) continue;
            for (const auto &entry : cols_to_rows[static_cast<size_t>(j)])
                subgrad[static_cast<size_t>(entry.row)] -= entry.val;
        }

        double norm_sq = 0.0;
        for (int i = 0; i < nrows; ++i)
            norm_sq += subgrad[static_cast<size_t>(i)] * subgrad[static_cast<size_t>(i)];

        if (norm_sq < 1e-12) break; // optimal (Lagrangian = LP bound)

        // Step size: t = lambda * (UB - L(u)) / ||s||^2
        const double step = lambda * (ub - lagrangian_bound) / norm_sq;
        if (step < 1e-12) break;

        // Update multipliers: u_i = max(0, u_i + t * s_i)
        for (int i = 0; i < nrows; ++i) {
            u[static_cast<size_t>(i)] = std::max(0.0,
                u[static_cast<size_t>(i)] + step * subgrad[static_cast<size_t>(i)]);
        }

        result.iterations = iter + 1;
    }

    result.best_bound = best_bound;
    result.best_multipliers = u;

    if (verbosity >= 2) {
        const double elapsed = std::chrono::duration<double>(Clock::now() - start).count();
        fprintf(stderr, "  Lagrangian relaxation: %d iters, bound=%.6g, heuristic=%.6g (%.1fs)\n",
                result.iterations, result.best_bound,
                std::isfinite(result.best_heuristic_obj) ? result.best_heuristic_obj : -1.0,
                elapsed);
    }

    return result;
}

LagrangianResult lagrangian_relaxation_at_node(
    const BaseRelaxationModel &base,
    double incumbent_obj,
    const std::vector<double> &init_multipliers,
    const std::vector<BranchDecision> &decisions,
    int max_iterations,
    double time_limit) {

    using Clock = std::chrono::steady_clock;
    const auto start = Clock::now();

    const int nrows = base.nrows;
    const int ncols = base.ncols;
    LagrangianResult result;

    if (nrows <= 0 || ncols <= 0) return result;

    // Build fixing map from decisions: +1 = fixed to 1, -1 = fixed to 0
    std::vector<int8_t> col_fix(static_cast<size_t>(ncols), 0);
    for (const auto &d : decisions) {
        if (d.var_index >= 0 && d.var_index < ncols)
            col_fix[static_cast<size_t>(d.var_index)] =
                static_cast<int8_t>(d.fix_value == 1 ? 1 : -1);
    }

    // Initialize multipliers from LP duals (clamped >= 0)
    std::vector<double> u(static_cast<size_t>(nrows), 0.0);
    if (static_cast<int>(init_multipliers.size()) >= nrows) {
        for (int i = 0; i < nrows; ++i)
            u[static_cast<size_t>(i)] = std::max(0.0, init_multipliers[static_cast<size_t>(i)]);
    }

    const auto &cols_to_rows = base.cols_to_rows;
    if (static_cast<int>(cols_to_rows.size()) < ncols) return result;

    double lambda = 2.0;
    int no_improve_count = 0;
    const int halve_interval = 15; // faster halving for shorter runs
    double best_bound = -std::numeric_limits<double>::infinity();
    double ub = incumbent_obj;

    std::vector<double> x(static_cast<size_t>(ncols));
    std::vector<double> rc(static_cast<size_t>(ncols));
    std::vector<double> subgrad(static_cast<size_t>(nrows));
    std::vector<double> coverage(static_cast<size_t>(nrows));

    for (int iter = 0; iter < max_iterations; ++iter) {
        const double elapsed = std::chrono::duration<double>(Clock::now() - start).count();
        if (elapsed > time_limit) break;

        // --- Lagrangian subproblem with fixings ---
        double lagrangian_bound = 0.0;
        for (int i = 0; i < nrows; ++i)
            lagrangian_bound += u[static_cast<size_t>(i)] * base.rhs[static_cast<size_t>(i)];

        for (int j = 0; j < ncols; ++j) {
            const auto sj = static_cast<size_t>(j);
            double r = base.obj[sj];
            for (const auto &entry : cols_to_rows[sj])
                r -= u[static_cast<size_t>(entry.row)] * entry.val;
            rc[sj] = r;

            if (col_fix[sj] == 1) {
                // Forced to 1
                x[sj] = 1.0;
                lagrangian_bound += r;
            } else if (col_fix[sj] == -1) {
                // Forced to 0
                x[sj] = 0.0;
            } else {
                // Free: standard Lagrangian decision
                if (r < 0.0) { x[sj] = 1.0; lagrangian_bound += r; }
                else { x[sj] = 0.0; }
            }
        }

        if (lagrangian_bound > best_bound) {
            best_bound = lagrangian_bound;
            no_improve_count = 0;
        } else {
            ++no_improve_count;
        }

        if (no_improve_count >= halve_interval) {
            lambda *= 0.5;
            no_improve_count = 0;
            if (lambda < 1e-6) break;
        }

        // --- Heuristic: seed from Lagrangian, repair, drop ---
        std::fill(coverage.begin(), coverage.end(), 0.0);
        std::vector<double> heur_sol(static_cast<size_t>(ncols), 0.0);

        // Include forced-to-1 and Lagrangian-selected free columns
        for (int j = 0; j < ncols; ++j) {
            const auto sj = static_cast<size_t>(j);
            if (col_fix[sj] == 1 || x[sj] > 0.5) {
                heur_sol[sj] = 1.0;
                for (const auto &entry : cols_to_rows[sj])
                    coverage[static_cast<size_t>(entry.row)] += entry.val;
            }
        }

        // Greedy repair: skip fixed-to-0 columns
        while (true) {
            int uncovered = -1;
            for (int i = 0; i < nrows; ++i) {
                if (coverage[static_cast<size_t>(i)] + 1e-8 < base.rhs[static_cast<size_t>(i)]) {
                    uncovered = i;
                    break;
                }
            }
            if (uncovered < 0) break;

            int best_col = -1;
            double best_score = -std::numeric_limits<double>::infinity();
            for (int j = 0; j < ncols; ++j) {
                const auto sj = static_cast<size_t>(j);
                if (heur_sol[sj] > 0.5 || col_fix[sj] == -1) continue;
                double uncov_gain = 0.0;
                for (const auto &entry : cols_to_rows[sj]) {
                    if (coverage[static_cast<size_t>(entry.row)] + 1e-8 <
                        base.rhs[static_cast<size_t>(entry.row)] && entry.val > 0.0)
                        uncov_gain += entry.val;
                }
                if (uncov_gain <= 0.0) continue;
                const double cost = std::max(1e-9, base.obj[sj]);
                const double score = uncov_gain / cost;
                if (score > best_score) {
                    best_score = score;
                    best_col = j;
                }
            }

            if (best_col < 0) break;
            heur_sol[static_cast<size_t>(best_col)] = 1.0;
            for (const auto &entry : cols_to_rows[static_cast<size_t>(best_col)])
                coverage[static_cast<size_t>(entry.row)] += entry.val;
        }

        // Drop redundant columns (most expensive first), skip forced-to-1
        std::vector<int> selected;
        for (int j = 0; j < ncols; ++j) {
            if (heur_sol[static_cast<size_t>(j)] > 0.5 && col_fix[static_cast<size_t>(j)] != 1)
                selected.push_back(j);
        }
        std::sort(selected.begin(), selected.end(),
                  [&](int a, int b) {
                      return base.obj[static_cast<size_t>(a)] > base.obj[static_cast<size_t>(b)];
                  });

        for (int col : selected) {
            bool can_remove = true;
            for (const auto &entry : cols_to_rows[static_cast<size_t>(col)]) {
                if (coverage[static_cast<size_t>(entry.row)] - entry.val + 1e-8 <
                    base.rhs[static_cast<size_t>(entry.row)]) {
                    can_remove = false;
                    break;
                }
            }
            if (can_remove) {
                heur_sol[static_cast<size_t>(col)] = 0.0;
                for (const auto &entry : cols_to_rows[static_cast<size_t>(col)])
                    coverage[static_cast<size_t>(entry.row)] -= entry.val;
            }
        }

        // Verify feasibility against CSR matrix
        bool feasible = true;
        for (int i = 0; i < nrows; ++i) {
            double row_sum = 0.0;
            for (int k = base.csr_offs[static_cast<size_t>(i)];
                 k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
                const int c = base.csr_inds[static_cast<size_t>(k)];
                if (c >= 0 && c < ncols && heur_sol[static_cast<size_t>(c)] > 0.5)
                    row_sum += base.csr_vals[static_cast<size_t>(k)];
            }
            if (row_sum + 1e-8 < base.rhs[static_cast<size_t>(i)]) {
                feasible = false;
                break;
            }
        }

        if (!feasible) continue;

        double heur_obj = 0.0;
        for (int j = 0; j < ncols; ++j) {
            if (heur_sol[static_cast<size_t>(j)] > 0.5)
                heur_obj += base.obj[static_cast<size_t>(j)];
        }

        if (heur_obj < result.best_heuristic_obj) {
            result.best_heuristic_obj = heur_obj;
            result.best_solution = heur_sol;
            if (heur_obj < ub)
                ub = heur_obj;
        }

        // --- Subgradient update ---
        std::fill(subgrad.begin(), subgrad.end(), 0.0);
        for (int i = 0; i < nrows; ++i)
            subgrad[static_cast<size_t>(i)] = base.rhs[static_cast<size_t>(i)];

        for (int j = 0; j < ncols; ++j) {
            if (x[static_cast<size_t>(j)] < 0.5) continue;
            for (const auto &entry : cols_to_rows[static_cast<size_t>(j)])
                subgrad[static_cast<size_t>(entry.row)] -= entry.val;
        }

        double norm_sq = 0.0;
        for (int i = 0; i < nrows; ++i)
            norm_sq += subgrad[static_cast<size_t>(i)] * subgrad[static_cast<size_t>(i)];

        if (norm_sq < 1e-12) break;

        const double step = lambda * (ub - lagrangian_bound) / norm_sq;
        if (step < 1e-12) break;

        for (int i = 0; i < nrows; ++i) {
            u[static_cast<size_t>(i)] = std::max(0.0,
                u[static_cast<size_t>(i)] + step * subgrad[static_cast<size_t>(i)]);
        }

        result.iterations = iter + 1;
    }

    result.best_bound = best_bound;
    result.best_multipliers = u;
    return result;
}

} // namespace scpsol
