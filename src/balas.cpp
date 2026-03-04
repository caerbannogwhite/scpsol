#include "balas.h"

#include <algorithm>
#include <cmath>

namespace scpsol {

BalasBranchResult balas_branch_generate(
    const std::vector<double> &primal_solution,
    const std::vector<double> &dual_solution,
    const std::vector<double> &reduced_costs,
    const BaseRelaxationModel &base,
    int ncols,
    int max_branches,
    double integrality_tol,
    bool force_balas) {
    BalasBranchResult result;
    result.use_br1 = false;

    std::vector<bool> is_fractional(static_cast<size_t>(ncols), false);
    for (int j = 0; j < ncols; ++j) {
        const double v = primal_solution[static_cast<size_t>(j)];
        const double nearest = std::floor(v + 0.5);
        if (std::fabs(v - nearest) > integrality_tol) {
            is_fractional[static_cast<size_t>(j)] = true;
        }
    }

    for (int i = 0; i < base.nrows; ++i) {
        if (dual_solution[static_cast<size_t>(i)] <= integrality_tol) continue;
        std::vector<int> Ri;
        for (int k = base.csr_offs[static_cast<size_t>(i)]; k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
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

    const int p = static_cast<int>(result.sets.size());
    if (p >= 2) {
        if (force_balas) {
            result.use_br1 = (p <= max_branches);
        } else {
            int total_vars = 0;
            int singleton_count = 0;
            for (const auto &s : result.sets) {
                total_vars += static_cast<int>(s.size());
                if (s.size() == 1) ++singleton_count;
            }
            result.use_br1 = (total_vars > static_cast<int>(p * std::log2(static_cast<double>(p)))) &&
                              (singleton_count <= 1) &&
                              (p <= max_branches);
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
