#include "diving.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

namespace scpsol {

// Collect fractional variable indices
static std::vector<int> get_fractional(
    const std::vector<double> &col_value, int ncols, double tol) {
    std::vector<int> frac;
    for (int j = 0; j < ncols; ++j) {
        const double v = col_value[static_cast<size_t>(j)];
        if (v > tol && v < 1.0 - tol)
            frac.push_back(j);
    }
    return frac;
}

// Greedy repair: given a partial integer solution, cover remaining rows
static bool greedy_repair(
    std::vector<double> &sol,
    const BaseRelaxationModel &base,
    double budget) {

    const int nrows = base.nrows;
    const int ncols = base.ncols;
    const auto &cols_to_rows = base.cols_to_rows;

    // Compute coverage
    std::vector<double> coverage(static_cast<size_t>(nrows), 0.0);
    double obj = 0.0;
    for (int j = 0; j < ncols; ++j) {
        if (sol[static_cast<size_t>(j)] > 0.5) {
            obj += base.obj[static_cast<size_t>(j)];
            for (const auto &e : cols_to_rows[static_cast<size_t>(j)])
                coverage[static_cast<size_t>(e.row)] += e.val;
        }
    }

    // Greedily add columns to cover uncovered rows
    while (true) {
        int uncov = -1;
        for (int i = 0; i < nrows; ++i) {
            if (coverage[static_cast<size_t>(i)] + 1e-8 < base.rhs[static_cast<size_t>(i)]) {
                uncov = i;
                break;
            }
        }
        if (uncov < 0) break; // all covered

        int best_col = -1;
        double best_score = -1.0;
        for (int j = 0; j < ncols; ++j) {
            if (sol[static_cast<size_t>(j)] > 0.5) continue;
            double gain = 0.0;
            for (const auto &e : cols_to_rows[static_cast<size_t>(j)]) {
                if (coverage[static_cast<size_t>(e.row)] + 1e-8 <
                    base.rhs[static_cast<size_t>(e.row)])
                    gain += 1.0;
            }
            if (gain <= 0.0) continue;
            double cost = std::max(1e-9, base.obj[static_cast<size_t>(j)]);
            double score = gain / cost;
            if (score > best_score) {
                best_score = score;
                best_col = j;
            }
        }
        if (best_col < 0) return false; // infeasible

        sol[static_cast<size_t>(best_col)] = 1.0;
        obj += base.obj[static_cast<size_t>(best_col)];
        if (obj >= budget) return false; // exceeded budget
        for (const auto &e : cols_to_rows[static_cast<size_t>(best_col)])
            coverage[static_cast<size_t>(e.row)] += e.val;
    }
    return true;
}

// Redundancy removal: drop expensive columns if coverage remains
static void drop_redundant(
    std::vector<double> &sol,
    const BaseRelaxationModel &base) {

    const int ncols = base.ncols;
    const int nrows = base.nrows;
    const auto &cols_to_rows = base.cols_to_rows;

    std::vector<double> coverage(static_cast<size_t>(nrows), 0.0);
    std::vector<int> selected;
    for (int j = 0; j < ncols; ++j) {
        if (sol[static_cast<size_t>(j)] > 0.5) {
            selected.push_back(j);
            for (const auto &e : cols_to_rows[static_cast<size_t>(j)])
                coverage[static_cast<size_t>(e.row)] += e.val;
        }
    }

    // Sort by descending cost
    std::sort(selected.begin(), selected.end(),
              [&](int a, int b) {
                  return base.obj[static_cast<size_t>(a)] > base.obj[static_cast<size_t>(b)];
              });

    for (int col : selected) {
        bool can_remove = true;
        for (const auto &e : cols_to_rows[static_cast<size_t>(col)]) {
            if (coverage[static_cast<size_t>(e.row)] - e.val + 1e-8 <
                base.rhs[static_cast<size_t>(e.row)]) {
                can_remove = false;
                break;
            }
        }
        if (can_remove) {
            sol[static_cast<size_t>(col)] = 0.0;
            for (const auto &e : cols_to_rows[static_cast<size_t>(col)])
                coverage[static_cast<size_t>(e.row)] -= e.val;
        }
    }
}

// Compute objective of a binary solution
static double compute_obj(const std::vector<double> &sol,
                           const std::vector<double> &obj, int ncols) {
    double v = 0.0;
    for (int j = 0; j < ncols; ++j)
        if (sol[static_cast<size_t>(j)] > 0.5)
            v += obj[static_cast<size_t>(j)];
    return v;
}

// Verify feasibility against CSR matrix
static bool verify_feasible(const std::vector<double> &sol,
                             const BaseRelaxationModel &base) {
    for (int i = 0; i < base.nrows; ++i) {
        double row_sum = 0.0;
        for (int k = base.csr_offs[static_cast<size_t>(i)];
             k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int c = base.csr_inds[static_cast<size_t>(k)];
            if (c >= 0 && c < base.ncols && sol[static_cast<size_t>(c)] > 0.5)
                row_sum += base.csr_vals[static_cast<size_t>(k)];
        }
        if (row_sum + 1e-8 < base.rhs[static_cast<size_t>(i)])
            return false;
    }
    return true;
}

enum class DiveStrategy { Fractional, Coefficient, Guided };

// Single dive attempt with a given strategy.
static DivingResult dive_single(
    LpSolver &lp,
    const BaseRelaxationModel &base,
    const BranchNodeState &branch_node,
    const LpSolution &start_sol,
    double best_obj,
    double integ_tol,
    int max_lp_solves,
    DiveStrategy strategy) {

    DivingResult result;
    const int ncols = base.ncols;

    // Build accumulated decisions: start with node's decisions
    std::vector<BranchDecision> decisions = branch_node.decisions;
    LpSolution sol = start_sol;

    for (int dive_step = 0; dive_step < max_lp_solves; ++dive_step) {
        auto frac = get_fractional(sol.col_value, ncols, integ_tol);
        if (frac.empty()) {
            // Integer solution found!
            result.solution.assign(static_cast<size_t>(ncols), 0.0);
            for (int j = 0; j < ncols; ++j)
                result.solution[static_cast<size_t>(j)] =
                    (sol.col_value[static_cast<size_t>(j)] > 0.5) ? 1.0 : 0.0;

            double obj = compute_obj(result.solution, base.obj, ncols);
            if (obj < best_obj - 1e-6) {
                result.found = true;
                result.objective = obj;
            }
            result.lp_solves = dive_step;
            return result;
        }

        // Pick variable and rounding direction
        int pick = -1;
        int round_to = 0;

        switch (strategy) {
        case DiveStrategy::Fractional: {
            // Pick variable closest to an integer boundary
            double best_frac_dist = 1.0;
            for (int j : frac) {
                double v = sol.col_value[static_cast<size_t>(j)];
                double dist = std::min(v, 1.0 - v);
                if (dist < best_frac_dist) {
                    best_frac_dist = dist;
                    pick = j;
                    round_to = (v >= 0.5) ? 1 : 0;
                }
            }
            break;
        }
        case DiveStrategy::Coefficient: {
            // Pick cheapest fractional variable, round up (add cheap coverage)
            double best_cost = std::numeric_limits<double>::infinity();
            for (int j : frac) {
                double c = base.obj[static_cast<size_t>(j)];
                if (c < best_cost) {
                    best_cost = c;
                    pick = j;
                }
            }
            round_to = 1; // bias toward including cheap columns
            break;
        }
        case DiveStrategy::Guided: {
            // Round toward LP value (most fractional = closest to 0.5)
            // Pick the variable with largest fractionality for maximum progress
            double best_frac = 0.0;
            for (int j : frac) {
                double v = sol.col_value[static_cast<size_t>(j)];
                double f = 0.5 - std::abs(v - 0.5); // max at v=0.5
                if (f > best_frac) {
                    best_frac = f;
                    pick = j;
                    round_to = (v >= 0.5) ? 1 : 0;
                }
            }
            break;
        }
        }

        if (pick < 0) break;

        // Add fixing decision
        decisions.push_back({pick, round_to});

        // Apply all decisions, solve LP, then restore
        lp.apply_decisions(decisions);
        lp.add_cuts(branch_node.cuts);
        lp.restore_basis();
        LpSolution new_sol = lp.solve();
        ++result.lp_solves;
        lp.save_basis(); // warm-start next dive step from this basis
        lp.restore_base_state();

        if (!new_sol.solved || !new_sol.optimal || new_sol.infeasible) {
            // Dive hit infeasibility — try repair on last good solution
            break;
        }

        if (new_sol.primal_obj >= best_obj - 1e-6) {
            // Exceeded cutoff — abort dive
            break;
        }

        sol = new_sol;
    }

    // Dive ended without finding an integer solution.
    // Try greedy repair on the last LP solution.
    result.solution.assign(static_cast<size_t>(ncols), 0.0);
    for (int j = 0; j < ncols; ++j) {
        const double v = sol.col_value[static_cast<size_t>(j)];
        result.solution[static_cast<size_t>(j)] = (v > 0.5) ? 1.0 : 0.0;
    }

    if (greedy_repair(result.solution, base, best_obj)) {
        drop_redundant(result.solution, base);
        if (verify_feasible(result.solution, base)) {
            double obj = compute_obj(result.solution, base.obj, ncols);
            if (obj < best_obj - 1e-6) {
                result.found = true;
                result.objective = obj;
            }
        }
    }

    return result;
}

DivingResult run_diving_heuristics(
    LpSolver &lp,
    const BaseRelaxationModel &base,
    const BranchNodeState &branch_node,
    const LpSolution &current_sol,
    double best_obj,
    double integ_tol,
    int max_lp_solves) {

    DivingResult best;
    const int per_strategy_budget = max_lp_solves;

    // Save basis state before diving
    lp.save_basis();

    DiveStrategy strategies[] = {
        DiveStrategy::Coefficient,
        DiveStrategy::Fractional,
        DiveStrategy::Guided,
    };
    const char *names[] = {"dive_coefficient", "dive_fractional", "dive_guided"};

    for (int s = 0; s < 3; ++s) {
        // Restore basis to pre-dive state for each strategy
        lp.restore_basis();

        auto dr = dive_single(lp, base, branch_node, current_sol,
                               best_obj, integ_tol, per_strategy_budget,
                               strategies[s]);

        if (dr.found && dr.objective < best.objective) {
            best = dr;
            best.strategy_name = names[s];
            best_obj = dr.objective; // tighten cutoff for subsequent strategies
        }
        best.lp_solves += dr.lp_solves;
    }

    // Restore basis to pre-dive state for BnB continuation
    lp.restore_basis();

    return best;
}

} // namespace scpsol
