#include "bnb.h"
#include "cuts.h"
#include "preprocessor.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>

namespace scpsol {

ModelReductionResult reduce_base_model(
    BaseRelaxationModel &base, double incumbent_bound, double tol) {
    ModelReductionResult result;
    result.old_to_new.assign(static_cast<size_t>(base.ncols), 0);

    std::vector<int> new_to_old;
    std::vector<int> new_active_to_original;
    int new_col = 0;
    for (int old_col = 0; old_col < base.ncols; ++old_col) {
        if (base.obj[static_cast<size_t>(old_col)] + tol >= incumbent_bound) {
            result.old_to_new[static_cast<size_t>(old_col)] = -1;
            ++result.columns_removed;
        } else {
            result.old_to_new[static_cast<size_t>(old_col)] = new_col;
            new_to_old.push_back(old_col);
            new_active_to_original.push_back(base.active_to_original[static_cast<size_t>(old_col)]);
            ++new_col;
        }
    }

    if (result.columns_removed == 0) return result;

    const int new_ncols = new_col;

    // Rebuild objective
    std::vector<double> new_obj(static_cast<size_t>(new_ncols));
    for (int j = 0; j < new_ncols; ++j)
        new_obj[static_cast<size_t>(j)] = base.obj[static_cast<size_t>(new_to_old[static_cast<size_t>(j)])];

    // Rebuild CSR
    std::vector<int> new_csr_inds;
    std::vector<int> new_csr_offs;
    std::vector<double> new_csr_vals;
    new_csr_offs.reserve(static_cast<size_t>(base.nrows) + 1);
    new_csr_offs.push_back(0);

    for (int i = 0; i < base.nrows; ++i) {
        for (int k = base.csr_offs[static_cast<size_t>(i)]; k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int old_col = base.csr_inds[static_cast<size_t>(k)];
            if (old_col >= 0 && old_col < base.ncols) {
                const int mapped = result.old_to_new[static_cast<size_t>(old_col)];
                if (mapped >= 0) {
                    new_csr_inds.push_back(mapped);
                    new_csr_vals.push_back(base.csr_vals[static_cast<size_t>(k)]);
                }
            }
        }
        new_csr_offs.push_back(static_cast<int>(new_csr_vals.size()));
    }

    base.ncols = new_ncols;
    base.nnz = static_cast<int>(new_csr_vals.size());
    base.csr_inds = std::move(new_csr_inds);
    base.csr_offs = std::move(new_csr_offs);
    base.csr_vals = std::move(new_csr_vals);
    base.obj = std::move(new_obj);
    base.active_to_original = std::move(new_active_to_original);

    // Remap base_cuts
    for (CutConstraint &cut : base.base_cuts) {
        remap_cut_constraint(cut, result.old_to_new);
    }

    return result;
}

ModelReductionResult reduce_base_model_budget_pruning(
    BaseRelaxationModel &base, double incumbent_bound, double tol,
    double preprocess_time_limit_sec) {
    ModelReductionResult result;
    result.old_to_new.assign(static_cast<size_t>(base.ncols), 0);

    if (base.nrows <= 0 || base.ncols <= 0 || !std::isfinite(incumbent_bound))
        return result;

    ColumnPreprocessContext ctx;
    ctx.nrows = base.nrows;
    ctx.ncols = base.ncols;
    ctx.costs.assign(base.obj.begin(), base.obj.begin() + base.ncols);
    ctx.active.assign(static_cast<size_t>(ctx.ncols), 1);
    ctx.incumbent_bound = incumbent_bound;
    if (preprocess_time_limit_sec > 0.0) {
        ctx.deadline = std::chrono::steady_clock::now() +
                       std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                           std::chrono::duration<double>(preprocess_time_limit_sec));
    }

    ctx.rows_by_column.assign(static_cast<size_t>(ctx.ncols), std::vector<int>());
    for (int i = 0; i < base.nrows; ++i) {
        for (int k = base.csr_offs[static_cast<size_t>(i)]; k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int col = base.csr_inds[static_cast<size_t>(k)];
            if (col >= 0 && col < base.ncols && base.csr_vals[static_cast<size_t>(k)] > tol) {
                ctx.rows_by_column[static_cast<size_t>(col)].push_back(i);
            }
        }
    }

    auto rules = make_preprocess_rules("incumbent_budget");
    int removed_by_rules = 0;
    for (const auto &rule : rules) {
        removed_by_rules += rule->apply(ctx, tol);
    }

    if (removed_by_rules <= 0) return result;

    std::vector<int> new_to_old;
    std::vector<int> new_active_to_original;
    int new_col = 0;
    for (int old_col = 0; old_col < base.ncols; ++old_col) {
        if (!ctx.active[static_cast<size_t>(old_col)]) {
            result.old_to_new[static_cast<size_t>(old_col)] = -1;
            ++result.columns_removed;
        } else {
            result.old_to_new[static_cast<size_t>(old_col)] = new_col;
            new_to_old.push_back(old_col);
            new_active_to_original.push_back(base.active_to_original[static_cast<size_t>(old_col)]);
            ++new_col;
        }
    }

    const int new_ncols = new_col;

    std::vector<double> new_obj(static_cast<size_t>(new_ncols));
    for (int j = 0; j < new_ncols; ++j)
        new_obj[static_cast<size_t>(j)] = base.obj[static_cast<size_t>(new_to_old[static_cast<size_t>(j)])];

    std::vector<int> new_csr_inds;
    std::vector<int> new_csr_offs;
    std::vector<double> new_csr_vals;
    new_csr_offs.reserve(static_cast<size_t>(base.nrows) + 1);
    new_csr_offs.push_back(0);

    for (int i = 0; i < base.nrows; ++i) {
        for (int k = base.csr_offs[static_cast<size_t>(i)]; k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int old_col = base.csr_inds[static_cast<size_t>(k)];
            if (old_col >= 0 && old_col < base.ncols) {
                const int mapped = result.old_to_new[static_cast<size_t>(old_col)];
                if (mapped >= 0) {
                    new_csr_inds.push_back(mapped);
                    new_csr_vals.push_back(base.csr_vals[static_cast<size_t>(k)]);
                }
            }
        }
        new_csr_offs.push_back(static_cast<int>(new_csr_vals.size()));
    }

    base.ncols = new_ncols;
    base.nnz = static_cast<int>(new_csr_vals.size());
    base.csr_inds = std::move(new_csr_inds);
    base.csr_offs = std::move(new_csr_offs);
    base.csr_vals = std::move(new_csr_vals);
    base.obj = std::move(new_obj);
    base.active_to_original = std::move(new_active_to_original);

    for (CutConstraint &cut : base.base_cuts) {
        remap_cut_constraint(cut, result.old_to_new);
    }

    return result;
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
