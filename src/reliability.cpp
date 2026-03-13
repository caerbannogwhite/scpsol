#include "reliability.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

namespace scpsol {

// ================================================================
// PseudocostState
// ================================================================

void PseudocostState::init(int ncols) {
    sum_down.assign(static_cast<size_t>(ncols), 0.0);
    sum_up.assign(static_cast<size_t>(ncols), 0.0);
    count_down.assign(static_cast<size_t>(ncols), 0);
    count_up.assign(static_cast<size_t>(ncols), 0);
}

void PseudocostState::update_down(int j, double frac, double gain) {
    if (frac < 1e-12) return;
    sum_down[static_cast<size_t>(j)] += gain / frac;
    ++count_down[static_cast<size_t>(j)];
}

void PseudocostState::update_up(int j, double frac, double gain) {
    const double one_minus_frac = 1.0 - frac;
    if (one_minus_frac < 1e-12) return;
    sum_up[static_cast<size_t>(j)] += gain / one_minus_frac;
    ++count_up[static_cast<size_t>(j)];
}

double PseudocostState::estimate_down(int j, double frac, double obj_j) const {
    const int cnt = count_down[static_cast<size_t>(j)];
    if (cnt > 0)
        return frac * (sum_down[static_cast<size_t>(j)] / cnt);
    // Fallback: cost-proportional estimate
    return frac * obj_j;
}

double PseudocostState::estimate_up(int j, double frac, double obj_j) const {
    const int cnt = count_up[static_cast<size_t>(j)];
    const double one_minus_frac = 1.0 - frac;
    if (cnt > 0)
        return one_minus_frac * (sum_up[static_cast<size_t>(j)] / cnt);
    return one_minus_frac * obj_j;
}

bool PseudocostState::is_reliable(int j, int eta_rel) const {
    return count_down[static_cast<size_t>(j)] >= eta_rel &&
           count_up[static_cast<size_t>(j)] >= eta_rel;
}

void PseudocostState::remap(const std::vector<int> &old_to_new, int new_ncols) {
    std::vector<double> new_sd(static_cast<size_t>(new_ncols), 0.0);
    std::vector<double> new_su(static_cast<size_t>(new_ncols), 0.0);
    std::vector<int> new_cd(static_cast<size_t>(new_ncols), 0);
    std::vector<int> new_cu(static_cast<size_t>(new_ncols), 0);
    for (int j = 0; j < static_cast<int>(old_to_new.size()); ++j) {
        const int nj = old_to_new[static_cast<size_t>(j)];
        if (nj >= 0 && nj < new_ncols) {
            new_sd[static_cast<size_t>(nj)] = sum_down[static_cast<size_t>(j)];
            new_su[static_cast<size_t>(nj)] = sum_up[static_cast<size_t>(j)];
            new_cd[static_cast<size_t>(nj)] = count_down[static_cast<size_t>(j)];
            new_cu[static_cast<size_t>(nj)] = count_up[static_cast<size_t>(j)];
        }
    }
    sum_down = std::move(new_sd);
    sum_up = std::move(new_su);
    count_down = std::move(new_cd);
    count_up = std::move(new_cu);
}

// ================================================================
// Adjacency and propagation
// ================================================================

ScpAdjacency build_adjacency(const BaseRelaxationModel &base) {
    ScpAdjacency adj;
    adj.nrows = base.nrows;
    adj.ncols = base.ncols;
    adj.cols_by_row.resize(static_cast<size_t>(base.nrows));
    adj.rows_by_col.resize(static_cast<size_t>(base.ncols));
    for (int i = 0; i < base.nrows; ++i) {
        for (int k = base.csr_offs[static_cast<size_t>(i)];
             k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int c = base.csr_inds[static_cast<size_t>(k)];
            if (c >= 0 && c < base.ncols) {
                adj.cols_by_row[static_cast<size_t>(i)].push_back(c);
                adj.rows_by_col[static_cast<size_t>(c)].push_back(i);
            }
        }
    }
    return adj;
}

PropagationResult propagate_fixings(
    const ScpAdjacency &adj,
    const std::vector<BranchDecision> &current_decisions,
    int probe_var, int probe_value) {

    PropagationResult result;
    const int ncols = adj.ncols;
    const int nrows = adj.nrows;

    // Build fixed state: 0 = free, 1 = fixed-to-1, -1 = fixed-to-0
    std::vector<signed char> col_state(static_cast<size_t>(ncols), 0);
    for (const auto &d : current_decisions) {
        if (d.var_index >= 0 && d.var_index < ncols)
            col_state[static_cast<size_t>(d.var_index)] =
                static_cast<signed char>(d.fix_value == 1 ? 1 : -1);
    }
    if (probe_var >= 0 && probe_var < ncols)
        col_state[static_cast<size_t>(probe_var)] =
            static_cast<signed char>(probe_value == 1 ? 1 : -1);

    // Track which rows are already covered (by a column fixed to 1)
    std::vector<char> row_covered(static_cast<size_t>(nrows), 0);
    for (int j = 0; j < ncols; ++j) {
        if (col_state[static_cast<size_t>(j)] == 1) {
            for (int r : adj.rows_by_col[static_cast<size_t>(j)])
                row_covered[static_cast<size_t>(r)] = 1;
        }
    }

    // Essential column cascade
    bool changed = true;
    while (changed) {
        changed = false;
        for (int i = 0; i < nrows; ++i) {
            if (row_covered[static_cast<size_t>(i)]) continue;
            int free_count = 0, last_free = -1;
            for (int c : adj.cols_by_row[static_cast<size_t>(i)]) {
                const auto st = col_state[static_cast<size_t>(c)];
                if (st == 1) {
                    // Row already covered — mark and skip
                    row_covered[static_cast<size_t>(i)] = 1;
                    free_count = -1;
                    break;
                }
                if (st == 0) { // free
                    ++free_count;
                    last_free = c;
                    if (free_count > 1) break; // no need to count further
                }
            }
            if (free_count == -1) continue; // row was covered
            if (free_count == 0) {
                result.infeasible = true;
                return result;
            }
            if (free_count == 1) {
                // Force last_free to 1
                col_state[static_cast<size_t>(last_free)] = 1;
                result.implied.push_back({last_free, 1});
                for (int r : adj.rows_by_col[static_cast<size_t>(last_free)])
                    row_covered[static_cast<size_t>(r)] = 1;
                changed = true;
            }
        }
    }
    return result;
}

// ================================================================
// Reliability branching selection
// ================================================================

int reliability_branch_select(
    PseudocostState &pc_state,
    LpSolver &lp,
    const BaseRelaxationModel &base,
    const BranchNodeState &branch_node,
    const LpSolution &current_sol,
    const std::vector<int> &fractional_candidates,
    double best_obj,
    int eta_rel,
    int max_sb_candidates,
    double /*integ_tol*/,
    int verbosity,
    int &lp_solves_out) {

    if (fractional_candidates.empty()) return -1;

    const double parent_dual = current_sol.dual_obj;
    const double eps = 1e-6;
    lp_solves_out = 0;

    // Build adjacency once for propagation
    ScpAdjacency adj = build_adjacency(base);

    // Partition candidates: unreliable first (need strong branching), then reliable.
    // Among unreliable, prioritize those closest to 0.5 (most fractional).
    std::vector<int> candidates = fractional_candidates;
    std::sort(candidates.begin(), candidates.end(), [&](int a, int b) {
        const bool ra = pc_state.is_reliable(a, eta_rel);
        const bool rb = pc_state.is_reliable(b, eta_rel);
        if (ra != rb) return !ra; // unreliable first
        const double fa = std::fabs(current_sol.col_value[static_cast<size_t>(a)] - 0.5);
        const double fb = std::fabs(current_sol.col_value[static_cast<size_t>(b)] - 0.5);
        return fa < fb; // most fractional first
    });

    int best_var = -1;
    double best_score = -1.0;
    int sb_count = 0;

    for (int j : candidates) {
        const double xj = current_sol.col_value[static_cast<size_t>(j)];
        const double frac = xj;       // for binary: frac in (0,1)
        const double obj_j = base.obj[static_cast<size_t>(j)];
        double d_down, d_up;

        if (pc_state.is_reliable(j, eta_rel) || sb_count >= max_sb_candidates) {
            // Use pseudocost estimates
            d_down = pc_state.estimate_down(j, frac, obj_j);
            d_up = pc_state.estimate_up(j, frac, obj_j);
        } else {
            // Strong branch: probe both directions
            ++sb_count;

            // --- Probe x_j = 0 ---
            auto prop0 = propagate_fixings(adj, branch_node.decisions, j, 0);
            if (prop0.infeasible) {
                d_down = 1e20;
            } else {
                // Build combined decisions: parent + probe + implied
                std::vector<BranchDecision> probe_decs = branch_node.decisions;
                probe_decs.push_back({j, 0});
                for (const auto &imp : prop0.implied)
                    probe_decs.push_back(imp);

                lp.apply_decisions(probe_decs);
                lp.add_cuts(branch_node.cuts);
                lp.restore_basis();
                LpSolution probe_sol = lp.solve();
                ++lp_solves_out;
                lp.restore_base_state();

                if (!probe_sol.solved || probe_sol.infeasible) {
                    d_down = 1e20;
                } else {
                    d_down = std::max(0.0, probe_sol.dual_obj - parent_dual);
                }
            }

            // --- Probe x_j = 1 ---
            auto prop1 = propagate_fixings(adj, branch_node.decisions, j, 1);
            if (prop1.infeasible) {
                d_up = 1e20;
            } else {
                std::vector<BranchDecision> probe_decs = branch_node.decisions;
                probe_decs.push_back({j, 1});
                for (const auto &imp : prop1.implied)
                    probe_decs.push_back(imp);

                lp.apply_decisions(probe_decs);
                lp.add_cuts(branch_node.cuts);
                lp.restore_basis();
                LpSolution probe_sol = lp.solve();
                ++lp_solves_out;
                lp.restore_base_state();

                if (!probe_sol.solved || probe_sol.infeasible) {
                    d_up = 1e20;
                } else {
                    d_up = std::max(0.0, probe_sol.dual_obj - parent_dual);
                }
            }

            // Both infeasible => current node infeasible
            if (d_down >= 1e19 && d_up >= 1e19)
                return -1;

            // Update pseudocosts (only for finite gains)
            if (d_down < 1e19)
                pc_state.update_down(j, frac, d_down);
            if (d_up < 1e19)
                pc_state.update_up(j, frac, d_up);

            // If one direction infeasible, the variable is forced.
            // We still score it normally — the product score will be huge,
            // making it the top candidate. The child creation in solver.cpp
            // will naturally prune the infeasible branch.
        }

        // Product score: max(eps, d_down) * max(eps, d_up)
        const double score = std::max(eps, d_down) * std::max(eps, d_up);
        if (score > best_score) {
            best_score = score;
            best_var = j;
        }
    }

    if (verbosity >= 4 && sb_count > 0)
        fprintf(stderr, "    Reliability branching: %d SB probes, %d LP solves, "
                        "best_var=%d score=%.6g\n",
                sb_count, lp_solves_out, best_var, best_score);

    return best_var;
}

} // namespace scpsol
