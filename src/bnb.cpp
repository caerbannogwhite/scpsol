#include "bnb.h"
#include "compact.h"
#include "cuts.h"
#include "preprocessor.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>

namespace scpsol {

namespace {

// Drop the columns with keep[j] == 0 from the base model (matrix, costs,
// column map and cuts) and return the old->new map.
ModelReductionResult remove_base_columns(BaseRelaxationModel &base,
                                         const std::vector<char> &keep) {
    ModelReductionResult result;
    const int new_ncols = build_column_map(keep, result.old_to_new);
    result.columns_removed = base.ncols - new_ncols;
    if (result.columns_removed <= 0) return result;

    compact_csr_columns(base.nrows, result.old_to_new, base.csr_inds, base.csr_offs, base.csr_vals);
    compact_column_vector(result.old_to_new, new_ncols, base.obj);
    compact_column_vector(result.old_to_new, new_ncols, base.active_to_original);
    base.ncols = new_ncols;
    base.nnz = static_cast<int>(base.csr_vals.size());

    for (CutConstraint &cut : base.base_cuts) {
        remap_cut_constraint(cut, result.old_to_new);
    }
    return result;
}

} // anonymous namespace

ModelReductionResult reduce_base_model(
    BaseRelaxationModel &base, double incumbent_bound, double tol) {
    std::vector<char> keep(static_cast<size_t>(base.ncols), 1);
    for (int j = 0; j < base.ncols; ++j) {
        if (base.obj[static_cast<size_t>(j)] + tol >= incumbent_bound)
            keep[static_cast<size_t>(j)] = 0;
    }
    return remove_base_columns(base, keep);
}

ModelReductionResult reduce_base_model_budget_pruning(
    BaseRelaxationModel &base, double incumbent_bound, double tol,
    double preprocess_time_limit_sec) {
    if (base.nrows <= 0 || base.ncols <= 0 || !std::isfinite(incumbent_bound))
        return ModelReductionResult{};

    ColumnPreprocessContext ctx = make_column_context(
        base.nrows, base.ncols, base.csr_inds, base.csr_offs, base.csr_vals, base.obj,
        tol, preprocess_time_limit_sec, incumbent_bound);

    auto rules = make_preprocess_rules("incumbent_budget");
    int removed_by_rules = 0;
    for (const auto &rule : rules) {
        removed_by_rules += rule->apply(ctx, tol);
    }
    if (removed_by_rules <= 0) return ModelReductionResult{};

    return remove_base_columns(base, ctx.active);
}

bool remap_branch_node(BranchNodeState &node, const std::vector<int> &old_to_new) {
    std::vector<BranchDecision> new_decisions;
    new_decisions.reserve(node.decisions.size());
    for (const BranchDecision &d : node.decisions) {
        if (d.var_index < 0 || d.var_index >= static_cast<int>(old_to_new.size())) continue;
        const int new_var = old_to_new[static_cast<size_t>(d.var_index)];
        if (new_var < 0) {
            if (d.fix_value == 1) return false;
            continue;
        }
        new_decisions.push_back({new_var, d.fix_value});
    }
    node.decisions = std::move(new_decisions);

    for (CutConstraint &cut : node.cuts) {
        remap_cut_constraint(cut, old_to_new);
    }
    return true;
}

bool append_decision_if_consistent(
    const BranchNodeState &parent, int var, int value, BranchNodeState *child) {
    *child = parent;
    for (const BranchDecision &d : parent.decisions) {
        if (d.var_index == var) {
            return d.fix_value == value;
        }
    }
    child->decisions.push_back({var, value});
    child->depth = parent.depth + 1;
    return true;
}

bool is_binary_integral(const std::vector<double> &x, int ncols, double tol) {
    for (int j = 0; j < ncols; ++j) {
        const double v = x[static_cast<size_t>(j)];
        const double nearest = floor(v + 0.5);
        if (fabs(v - nearest) > tol) return false;
        if (nearest < -tol || nearest > 1.0 + tol) return false;
    }
    return true;
}

std::vector<int> collect_fractional_candidates(
    const std::vector<double> &x, int ncols, double tol) {
    std::vector<int> candidates;
    candidates.reserve(static_cast<size_t>(ncols));
    for (int j = 0; j < ncols; ++j) {
        const double v = x[static_cast<size_t>(j)];
        const double nearest = floor(v + 0.5);
        if (fabs(v - nearest) > tol || nearest < -tol || nearest > 1.0 + tol) {
            candidates.push_back(j);
        }
    }
    return candidates;
}

bool is_node_provably_infeasible(
    const BranchNodeState &node, const BaseRelaxationModel &base) {
    std::vector<int> zero_fixed;
    for (const BranchDecision &d : node.decisions) {
        if (d.fix_value == 0) zero_fixed.push_back(d.var_index);
    }
    if (zero_fixed.empty()) return false;
    std::sort(zero_fixed.begin(), zero_fixed.end());

    // Check node-level cuts
    for (const CutConstraint &cut : node.cuts) {
        if (cut.rhs <= 0.0) continue;
        bool has_unfixed = false;
        for (const int j : cut.indices) {
            if (!std::binary_search(zero_fixed.begin(), zero_fixed.end(), j)) {
                has_unfixed = true;
                break;
            }
        }
        if (!has_unfixed) return true;
    }

    // Check base model cover rows
    for (int i = 0; i < base.nrows; ++i) {
        if (base.rhs[static_cast<size_t>(i)] <= 0.0) continue;
        bool has_unfixed = false;
        for (int k = base.csr_offs[static_cast<size_t>(i)]; k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int col = base.csr_inds[static_cast<size_t>(k)];
            if (col < 0 || col >= base.ncols) continue;
            if (base.csr_vals[static_cast<size_t>(k)] <= 0.0) continue;
            if (!std::binary_search(zero_fixed.begin(), zero_fixed.end(), col)) {
                has_unfixed = true;
                break;
            }
        }
        if (!has_unfixed) return true;
    }

    return false;
}

bool has_integer_objective(const std::vector<double> &obj, int ncols, double tol) {
    for (int j = 0; j < ncols; ++j) {
        const double v = obj[static_cast<size_t>(j)];
        if (fabs(v - floor(v + 0.5)) > tol) return false;
    }
    return true;
}

double tighten_dual_bound(double bound, double tol) {
    if (!std::isfinite(bound)) return bound;
    return ceil(bound - tol);
}

double compute_mip_gap(double incumbent, double dual_bound) {
    if (!std::isfinite(incumbent) || !std::isfinite(dual_bound))
        return std::numeric_limits<double>::infinity();
    if (dual_bound > incumbent)
        return std::numeric_limits<double>::infinity();
    return (incumbent - dual_bound) / std::max(1.0, fabs(incumbent));
}

} // namespace scpsol
