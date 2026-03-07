#include "balas.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace scpsol {

// ============================================================
// Shared helpers for Algorithms 2 and 3
// ============================================================

// Transpose of CSR: column j -> list of rows that cover j
static std::vector<std::vector<int>> build_col_to_rows(
    const BaseRelaxationModel &base, int ncols, double tol) {
    std::vector<std::vector<int>> col_to_rows(static_cast<size_t>(ncols));
    for (int i = 0; i < base.nrows; ++i) {
        for (int k = base.csr_offs[static_cast<size_t>(i)];
             k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int col = base.csr_inds[static_cast<size_t>(k)];
            if (col >= 0 && col < ncols && base.csr_vals[static_cast<size_t>(k)] > tol)
                col_to_rows[static_cast<size_t>(col)].push_back(i);
        }
    }
    return col_to_rows;
}

// T(x) = rows where constraint is tight (A^i * x ~ rhs_i)
static std::vector<bool> compute_tight_rows(
    const std::vector<double> &x,
    const BaseRelaxationModel &base, int ncols, double tol) {
    std::vector<bool> tight(static_cast<size_t>(base.nrows), false);
    for (int i = 0; i < base.nrows; ++i) {
        double lhs = 0.0;
        for (int k = base.csr_offs[static_cast<size_t>(i)];
             k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int col = base.csr_inds[static_cast<size_t>(k)];
            if (col >= 0 && col < ncols)
                lhs += base.csr_vals[static_cast<size_t>(k)] * x[static_cast<size_t>(col)];
        }
        if (std::fabs(lhs - base.rhs[static_cast<size_t>(i)]) <= tol)
            tight[static_cast<size_t>(i)] = true;
    }
    return tight;
}

// y = sum_i u_i * rhs_i (dual objective value)
static double compute_dual_objective(
    const std::vector<double> &dual,
    const BaseRelaxationModel &base) {
    double y = 0.0;
    for (int i = 0; i < base.nrows; ++i) {
        y += dual[static_cast<size_t>(i)] * base.rhs[static_cast<size_t>(i)];
    }
    return y;
}

// ============================================================
// Algorithm 2: BCG (Balas Cut Generation)
// ============================================================

CutConstraint balas_cut_generate(
    const std::vector<double> &primal,
    const std::vector<double> &dual,
    const std::vector<double> &reduced_costs,
    double incumbent_bound,
    const BaseRelaxationModel &base,
    int ncols, double tol) {
    CutConstraint empty_cut;

    if (!std::isfinite(incumbent_bound)) return empty_cut;

    // Check S first (cheap) before building expensive data structures
    std::vector<bool> in_S(static_cast<size_t>(ncols), false);
    std::vector<int> S;
    for (int j = 0; j < ncols; ++j) {
        if (primal[static_cast<size_t>(j)] > tol &&
            reduced_costs[static_cast<size_t>(j)] > tol) {
            in_S[static_cast<size_t>(j)] = true;
            S.push_back(j);
        }
    }
    if (S.empty()) return empty_cut;

    const auto col_to_rows = build_col_to_rows(base, ncols, tol);
    const auto tight = compute_tight_rows(primal, base, ncols, tol);
    const double z_u = incumbent_bound;

    // W = columns for the cut, s = working reduced costs
    std::vector<bool> in_W(static_cast<size_t>(ncols), false);
    std::vector<double> s(reduced_costs.begin(), reduced_costs.end());
    double y = compute_dual_objective(dual, base);

    bool success = false;

    while (!S.empty()) {
        const double gap = z_u - y;
        if (gap <= tol) { success = true; break; }

        // v_t: smallest s_j in S that alone reaches gap; if none, use max s_j
        double v_t = -1.0;
        double max_sj = -std::numeric_limits<double>::infinity();
        double min_reaching = std::numeric_limits<double>::infinity();
        for (const int j : S) {
            const double sj = s[static_cast<size_t>(j)];
            if (sj > max_sj) max_sj = sj;
            if (sj >= gap - tol && sj < min_reaching) min_reaching = sj;
        }
        v_t = (min_reaching < std::numeric_limits<double>::infinity()) ? min_reaching : max_sj;
        if (v_t <= tol) break;

        // J = {j in S | s_j ~ v_t}, Q = {j in 0..ncols | s_j >= v_t}
        std::vector<int> J;
        std::vector<bool> in_Q(static_cast<size_t>(ncols), false);
        for (const int j : S) {
            if (std::fabs(s[static_cast<size_t>(j)] - v_t) <= tol)
                J.push_back(j);
        }
        for (int j = 0; j < ncols; ++j) {
            if (s[static_cast<size_t>(j)] >= v_t - tol)
                in_Q[static_cast<size_t>(j)] = true;
        }

        // M_J = rows covered by columns in J
        std::vector<bool> in_MJ(static_cast<size_t>(base.nrows), false);
        for (const int j : J) {
            for (const int row : col_to_rows[static_cast<size_t>(j)])
                in_MJ[static_cast<size_t>(row)] = true;
        }

        // Build set of columns in J per row for fast lookup
        // N_i = columns covering row i (from CSR)
        // i* = row in T(x) & M_J minimizing |(N_i \ Q) \ W|
        int best_row = -1;
        int best_cost = std::numeric_limits<int>::max();
        for (int i = 0; i < base.nrows; ++i) {
            if (!tight[static_cast<size_t>(i)] || !in_MJ[static_cast<size_t>(i)]) continue;
            int cost = 0;
            for (int k = base.csr_offs[static_cast<size_t>(i)];
                 k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
                const int col = base.csr_inds[static_cast<size_t>(k)];
                if (col < 0 || col >= ncols) continue;
                if (base.csr_vals[static_cast<size_t>(k)] <= tol) continue;
                if (!in_Q[static_cast<size_t>(col)] && !in_W[static_cast<size_t>(col)])
                    ++cost;
            }
            if (cost < best_cost) {
                best_cost = cost;
                best_row = i;
            }
        }
        // Fallback: any row in M_J
        if (best_row < 0) {
            for (int i = 0; i < base.nrows; ++i) {
                if (in_MJ[static_cast<size_t>(i)]) { best_row = i; break; }
            }
        }
        if (best_row < 0) break;

        // j* = a column in J that covers best_row
        int best_col = -1;
        {
            // Build set of columns in J for fast lookup
            std::vector<bool> in_J(static_cast<size_t>(ncols), false);
            for (const int j : J) in_J[static_cast<size_t>(j)] = true;

            for (int k = base.csr_offs[static_cast<size_t>(best_row)];
                 k < base.csr_offs[static_cast<size_t>(best_row) + 1]; ++k) {
                const int col = base.csr_inds[static_cast<size_t>(k)];
                if (col >= 0 && col < ncols && in_J[static_cast<size_t>(col)] &&
                    base.csr_vals[static_cast<size_t>(k)] > tol) {
                    best_col = col;
                    break;
                }
            }
        }
        if (best_col < 0) break;

        // W = W union (N_{best_row} \ Q)
        for (int k = base.csr_offs[static_cast<size_t>(best_row)];
             k < base.csr_offs[static_cast<size_t>(best_row) + 1]; ++k) {
            const int col = base.csr_inds[static_cast<size_t>(k)];
            if (col >= 0 && col < ncols && base.csr_vals[static_cast<size_t>(k)] > tol &&
                !in_Q[static_cast<size_t>(col)])
                in_W[static_cast<size_t>(col)] = true;
        }

        // y += s[j*]
        y += s[static_cast<size_t>(best_col)];

        if (y >= z_u - tol) { success = true; break; }

        // Remove j* from S
        in_S[static_cast<size_t>(best_col)] = false;
        S.erase(std::remove(S.begin(), S.end(), best_col), S.end());

        // Update: s_j -= s[j*] for j in N_{best_row} & Q
        const double s_star = s[static_cast<size_t>(best_col)];
        for (int k = base.csr_offs[static_cast<size_t>(best_row)];
             k < base.csr_offs[static_cast<size_t>(best_row) + 1]; ++k) {
            const int col = base.csr_inds[static_cast<size_t>(k)];
            if (col >= 0 && col < ncols && base.csr_vals[static_cast<size_t>(k)] > tol &&
                in_Q[static_cast<size_t>(col)])
                s[static_cast<size_t>(col)] -= s_star;
        }
    }

    if (!success) return empty_cut;

    // Build cut: sum_{j in W} x_j >= 1
    CutConstraint cut;
    cut.type = "balas_bcg";
    cut.rhs = 1.0;
    for (int j = 0; j < ncols; ++j) {
        if (in_W[static_cast<size_t>(j)]) {
            cut.indices.push_back(j);
            cut.values.push_back(1.0);
        }
    }
    if (cut.indices.empty()) return empty_cut;
    return cut;
}

// ============================================================
// Algorithm 3: BBG (Balas Branch Generation)
// ============================================================

BalasBranchResult balas_branch_generate(
    const std::vector<double> &primal_solution,
    const std::vector<double> &dual_solution,
    const std::vector<double> &reduced_costs,
    double incumbent_bound,
    const BaseRelaxationModel &base,
    int ncols,
    int max_branches,
    double integrality_tol,
    bool force_balas) {
    BalasBranchResult result;
    result.use_br1 = false;

    // Try Algorithm 3 (iterative conditional-bound procedure)
    if (std::isfinite(incumbent_bound)) {
        // Check S first (cheap) before building expensive data structures
        std::vector<bool> in_S(static_cast<size_t>(ncols), false);
        std::vector<int> S;
        for (int j = 0; j < ncols; ++j) {
            if (primal_solution[static_cast<size_t>(j)] > integrality_tol &&
                reduced_costs[static_cast<size_t>(j)] > integrality_tol) {
                in_S[static_cast<size_t>(j)] = true;
                S.push_back(j);
            }
        }

        if (!S.empty()) {
            const auto col_to_rows = build_col_to_rows(base, ncols, integrality_tol);
            const auto tight = compute_tight_rows(primal_solution, base, ncols, integrality_tol);
            const double z_u = incumbent_bound;
            std::vector<double> s(reduced_costs.begin(), reduced_costs.end());
            double y = compute_dual_objective(dual_solution, base);

            while (!S.empty()) {
                const double gap = z_u - y;
                if (gap <= integrality_tol) break;

                // v_t computation
                double v_t = -1.0;
                double max_sj = -std::numeric_limits<double>::infinity();
                double min_reaching = std::numeric_limits<double>::infinity();
                for (const int j : S) {
                    const double sj = s[static_cast<size_t>(j)];
                    if (sj > max_sj) max_sj = sj;
                    if (sj >= gap - integrality_tol && sj < min_reaching) min_reaching = sj;
                }
                v_t = (min_reaching < std::numeric_limits<double>::infinity()) ? min_reaching : max_sj;
                if (v_t <= integrality_tol) break;

                // J = {j in S | s_j ~ v_t}
                std::vector<int> J;
                std::vector<bool> in_Q(static_cast<size_t>(ncols), false);
                for (const int j : S) {
                    if (std::fabs(s[static_cast<size_t>(j)] - v_t) <= integrality_tol)
                        J.push_back(j);
                }
                // Q = {j in 0..ncols | s_j >= v_t}
                for (int j = 0; j < ncols; ++j) {
                    if (s[static_cast<size_t>(j)] >= v_t - integrality_tol)
                        in_Q[static_cast<size_t>(j)] = true;
                }

                // M_J = rows covered by columns in J
                std::vector<bool> in_MJ(static_cast<size_t>(base.nrows), false);
                for (const int j : J) {
                    for (const int row : col_to_rows[static_cast<size_t>(j)])
                        in_MJ[static_cast<size_t>(row)] = true;
                }

                // i* = row in T(x) & M_J minimizing |N_i \ Q|
                int best_row = -1;
                int best_cost = std::numeric_limits<int>::max();
                for (int i = 0; i < base.nrows; ++i) {
                    if (!tight[static_cast<size_t>(i)] || !in_MJ[static_cast<size_t>(i)]) continue;
                    int cost = 0;
                    for (int k = base.csr_offs[static_cast<size_t>(i)];
                         k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
                        const int col = base.csr_inds[static_cast<size_t>(k)];
                        if (col < 0 || col >= ncols) continue;
                        if (base.csr_vals[static_cast<size_t>(k)] <= integrality_tol) continue;
                        if (!in_Q[static_cast<size_t>(col)]) ++cost;
                    }
                    if (cost < best_cost) {
                        best_cost = cost;
                        best_row = i;
                    }
                }
                if (best_row < 0) {
                    for (int i = 0; i < base.nrows; ++i) {
                        if (in_MJ[static_cast<size_t>(i)]) { best_row = i; break; }
                    }
                }
                if (best_row < 0) break;

                // j* = a column in J that covers best_row
                int best_col = -1;
                {
                    std::vector<bool> in_J(static_cast<size_t>(ncols), false);
                    for (const int j : J) in_J[static_cast<size_t>(j)] = true;
                    for (int k = base.csr_offs[static_cast<size_t>(best_row)];
                         k < base.csr_offs[static_cast<size_t>(best_row) + 1]; ++k) {
                        const int col = base.csr_inds[static_cast<size_t>(k)];
                        if (col >= 0 && col < ncols && in_J[static_cast<size_t>(col)] &&
                            base.csr_vals[static_cast<size_t>(k)] > integrality_tol) {
                            best_col = col;
                            break;
                        }
                    }
                }
                if (best_col < 0) break;

                // Collect N_{i*} ∩ Q as branch set R_k
                std::vector<int> R_k;
                for (int k = base.csr_offs[static_cast<size_t>(best_row)];
                     k < base.csr_offs[static_cast<size_t>(best_row) + 1]; ++k) {
                    const int col = base.csr_inds[static_cast<size_t>(k)];
                    if (col >= 0 && col < ncols &&
                        base.csr_vals[static_cast<size_t>(k)] > integrality_tol &&
                        in_Q[static_cast<size_t>(col)])
                        R_k.push_back(col);
                }
                if (!R_k.empty()) {
                    std::sort(R_k.begin(), R_k.end());
                    result.sets.push_back(std::move(R_k));
                }

                // y += s[j*]
                y += s[static_cast<size_t>(best_col)];

                if (y >= z_u - integrality_tol) break;

                // Remove j* from S
                in_S[static_cast<size_t>(best_col)] = false;
                S.erase(std::remove(S.begin(), S.end(), best_col), S.end());

                // Update: s_j -= s[j*] for j in N_{best_row} & Q
                const double s_star = s[static_cast<size_t>(best_col)];
                for (int k = base.csr_offs[static_cast<size_t>(best_row)];
                     k < base.csr_offs[static_cast<size_t>(best_row) + 1]; ++k) {
                    const int col = base.csr_inds[static_cast<size_t>(k)];
                    if (col >= 0 && col < ncols && base.csr_vals[static_cast<size_t>(k)] > integrality_tol &&
                        in_Q[static_cast<size_t>(col)])
                        s[static_cast<size_t>(col)] -= s_star;
                }
            }

            // Deduplicate
            std::sort(result.sets.begin(), result.sets.end());
            result.sets.erase(std::unique(result.sets.begin(), result.sets.end()), result.sets.end());
        }
    }

    // Fall back to per-row heuristic if iterative procedure produced < 2 sets
    if (static_cast<int>(result.sets.size()) < 2) {
        result.sets.clear();
        std::vector<bool> is_fractional(static_cast<size_t>(ncols), false);
        for (int j = 0; j < ncols; ++j) {
            const double v = primal_solution[static_cast<size_t>(j)];
            const double nearest = std::floor(v + 0.5);
            if (std::fabs(v - nearest) > integrality_tol)
                is_fractional[static_cast<size_t>(j)] = true;
        }
        for (int i = 0; i < base.nrows; ++i) {
            if (dual_solution[static_cast<size_t>(i)] <= integrality_tol) continue;
            std::vector<int> Ri;
            for (int k = base.csr_offs[static_cast<size_t>(i)];
                 k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
                const int col = base.csr_inds[static_cast<size_t>(k)];
                if (col < 0 || col >= ncols) continue;
                if (!is_fractional[static_cast<size_t>(col)]) continue;
                if (reduced_costs[static_cast<size_t>(col)] >= -integrality_tol) continue;
                if (base.csr_vals[static_cast<size_t>(k)] <= integrality_tol) continue;
                Ri.push_back(col);
            }
            if (!Ri.empty()) {
                std::sort(Ri.begin(), Ri.end());
                result.sets.push_back(std::move(Ri));
            }
        }
        // Deduplicate
        std::sort(result.sets.begin(), result.sets.end());
        result.sets.erase(std::unique(result.sets.begin(), result.sets.end()), result.sets.end());
    }

    // BR1 criteria (3.11-3.13) with tightened thresholds
    const int p = static_cast<int>(result.sets.size());
    if (p >= 3) {
        if (force_balas) {
            // Even forced, limit children to avoid frontier explosion
            result.use_br1 = (p <= std::min(max_branches, 8));
        } else {
            int total_vars = 0;
            int singleton_count = 0;
            for (const auto &s : result.sets) {
                total_vars += static_cast<int>(s.size());
                if (s.size() == 1) ++singleton_count;
            }
            // Require sufficient variable diversity and limit children
            const int effective_max = std::min(max_branches, 12);
            result.use_br1 = (total_vars > static_cast<int>(p * std::log2(static_cast<double>(p)))) &&
                              (singleton_count <= 1) &&
                              (p <= effective_max) &&
                              (total_vars >= 2 * p);
        }
    }

    return result;
}

std::vector<BranchNodeState> balas_br1_children(
    const BranchNodeState &parent,
    const std::vector<std::vector<int>> &branch_sets,
    double parent_dual_bound,
    double parent_dual_bound_raw) {
    const int p = static_cast<int>(branch_sets.size());
    std::vector<BranchNodeState> children;
    children.reserve(static_cast<size_t>(p));

    for (int i = 0; i < p; ++i) {
        BranchNodeState child = parent;
        bool consistent = true;
        for (const int j : branch_sets[static_cast<size_t>(i)]) {
            bool already_fixed = false;
            for (const BranchDecision &d : child.decisions) {
                if (d.var_index == j) {
                    if (d.fix_value != 0) consistent = false;
                    already_fixed = true;
                    break;
                }
            }
            if (!consistent) break;
            if (!already_fixed) {
                child.decisions.push_back({j, 0});
            }
        }
        if (!consistent) continue;

        // Collect zero-fixed vars for feasibility checks
        std::vector<int> zero_fixed;
        for (const BranchDecision &d : child.decisions) {
            if (d.fix_value == 0) zero_fixed.push_back(d.var_index);
        }
        std::sort(zero_fixed.begin(), zero_fixed.end());

        // Add cover cuts for prior sets and check feasibility
        bool cover_feasible = true;
        for (int k = 0; k < i; ++k) {
            bool has_unfixed = false;
            for (const int j : branch_sets[static_cast<size_t>(k)]) {
                if (!std::binary_search(zero_fixed.begin(), zero_fixed.end(), j)) {
                    has_unfixed = true;
                    break;
                }
            }
            if (!has_unfixed) { cover_feasible = false; break; }
        }
        if (!cover_feasible) continue;

        for (int k = 0; k < i; ++k) {
            CutConstraint cut;
            cut.type = "balas_cover";
            cut.rhs = 1.0;
            for (const int j : branch_sets[static_cast<size_t>(k)]) {
                cut.indices.push_back(j);
                cut.values.push_back(1.0);
            }
            child.cuts.push_back(std::move(cut));
        }

        child.depth = parent.depth + 1;
        child.parent_dual_bound = parent_dual_bound;
        child.parent_dual_bound_raw = parent_dual_bound_raw;
        children.push_back(std::move(child));
    }

    return children;
}

} // namespace scpsol
