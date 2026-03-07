#include "cuts.h"
#include "balas.h"

#include <algorithm>
#include <cmath>

namespace scpsol {

std::vector<double> compute_reduced_costs(
    const std::vector<double> &objective,
    const std::vector<double> &dual_solution,
    const BaseRelaxationModel &base,
    int ncols) {
    std::vector<double> s(static_cast<size_t>(ncols));
    for (int j = 0; j < ncols; ++j) {
        s[static_cast<size_t>(j)] = objective[static_cast<size_t>(j)];
    }
    for (int i = 0; i < base.nrows; ++i) {
        const double u = dual_solution[static_cast<size_t>(i)];
        if (u <= 0.0) continue;
        for (int k = base.csr_offs[static_cast<size_t>(i)]; k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int col = base.csr_inds[static_cast<size_t>(k)];
            if (col >= 0 && col < ncols) {
                s[static_cast<size_t>(col)] -= u * base.csr_vals[static_cast<size_t>(k)];
            }
        }
    }
    return s;
}

namespace {

class DualAggregatedCgSeparator : public ICutSeparator {
public:
    std::vector<CutConstraint> separate(
        const std::vector<double> &primal_solution,
        const std::vector<double> &dual_solution,
        const BaseRelaxationModel &base,
        int ncols,
        double tol,
        double /*incumbent_bound*/) const override {
        std::vector<CutConstraint> cuts;
        if (static_cast<int>(dual_solution.size()) < base.nrows ||
            static_cast<int>(primal_solution.size()) < ncols)
            return cuts;

        std::vector<double> agg_coeffs(static_cast<size_t>(ncols), 0.0);
        double rhs_sum = 0.0;
        for (int i = 0; i < base.nrows; ++i) {
            const double u = std::max(0.0, dual_solution[static_cast<size_t>(i)]);
            if (u < tol) continue;
            rhs_sum += u * base.rhs[static_cast<size_t>(i)];
            for (int k = base.csr_offs[static_cast<size_t>(i)]; k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
                const int col = base.csr_inds[static_cast<size_t>(k)];
                if (col >= 0 && col < ncols)
                    agg_coeffs[static_cast<size_t>(col)] += u * base.csr_vals[static_cast<size_t>(k)];
            }
        }

        const double cut_rhs = std::ceil(rhs_sum - tol);
        if (cut_rhs <= tol) return cuts;

        const double f0 = rhs_sum - std::floor(rhs_sum);
        if (f0 < tol || f0 > 1.0 - tol) return cuts;

        CutConstraint cut;
        cut.type = "cg_dual_aggregated";
        cut.rhs = cut_rhs;
        double lhs_val = 0.0;
        for (int j = 0; j < ncols; ++j) {
            const double rounded = std::ceil(agg_coeffs[static_cast<size_t>(j)] - tol);
            if (rounded > tol) {
                cut.indices.push_back(j);
                cut.values.push_back(rounded);
                lhs_val += rounded * primal_solution[static_cast<size_t>(j)];
            }
        }
        if (lhs_val < cut_rhs - tol && !cut.indices.empty()) {
            cuts.push_back(std::move(cut));
        }
        return cuts;
    }
};

class RowPairCgSeparator : public ICutSeparator {
public:
    std::vector<CutConstraint> separate(
        const std::vector<double> &primal_solution,
        const std::vector<double> &dual_solution,
        const BaseRelaxationModel &base,
        int ncols,
        double tol,
        double /*incumbent_bound*/) const override {
        std::vector<CutConstraint> cuts;
        if (static_cast<int>(dual_solution.size()) < base.nrows ||
            static_cast<int>(primal_solution.size()) < ncols)
            return cuts;

        std::vector<int> active_rows;
        active_rows.reserve(static_cast<size_t>(base.nrows));
        for (int i = 0; i < base.nrows; ++i) {
            if (dual_solution[static_cast<size_t>(i)] > tol)
                active_rows.push_back(i);
        }
        std::sort(active_rows.begin(), active_rows.end(),
                  [&](int a, int b) { return dual_solution[static_cast<size_t>(a)] > dual_solution[static_cast<size_t>(b)]; });

        const int max_active_rows = std::min(static_cast<int>(active_rows.size()), 40);
        std::vector<double> agg_coeffs(static_cast<size_t>(ncols), 0.0);

        for (int ri = 0; ri < max_active_rows; ++ri) {
            for (int rj = ri + 1; rj < max_active_rows; ++rj) {
                const int i1 = active_rows[static_cast<size_t>(ri)];
                const int i2 = active_rows[static_cast<size_t>(rj)];
                const double u1 = dual_solution[static_cast<size_t>(i1)];
                const double u2 = dual_solution[static_cast<size_t>(i2)];

                const double rhs_agg = u1 * base.rhs[static_cast<size_t>(i1)] +
                                       u2 * base.rhs[static_cast<size_t>(i2)];
                const double f0 = rhs_agg - std::floor(rhs_agg);
                if (f0 < tol || f0 > 1.0 - tol) continue;

                const double cut_rhs = std::ceil(rhs_agg - tol);
                if (cut_rhs <= tol) continue;

                std::fill(agg_coeffs.begin(), agg_coeffs.end(), 0.0);
                for (int k = base.csr_offs[static_cast<size_t>(i1)]; k < base.csr_offs[static_cast<size_t>(i1) + 1]; ++k) {
                    const int col = base.csr_inds[static_cast<size_t>(k)];
                    if (col >= 0 && col < ncols)
                        agg_coeffs[static_cast<size_t>(col)] += u1 * base.csr_vals[static_cast<size_t>(k)];
                }
                for (int k = base.csr_offs[static_cast<size_t>(i2)]; k < base.csr_offs[static_cast<size_t>(i2) + 1]; ++k) {
                    const int col = base.csr_inds[static_cast<size_t>(k)];
                    if (col >= 0 && col < ncols)
                        agg_coeffs[static_cast<size_t>(col)] += u2 * base.csr_vals[static_cast<size_t>(k)];
                }

                CutConstraint cut;
                cut.type = "cg_row_pair";
                cut.rhs = cut_rhs;
                double lhs_val = 0.0;
                for (int j = 0; j < ncols; ++j) {
                    const double rounded = std::ceil(agg_coeffs[static_cast<size_t>(j)] - tol);
                    if (rounded > tol) {
                        cut.indices.push_back(j);
                        cut.values.push_back(rounded);
                        lhs_val += rounded * primal_solution[static_cast<size_t>(j)];
                    }
                }
                if (lhs_val < cut_rhs - tol && !cut.indices.empty()) {
                    cuts.push_back(std::move(cut));
                    if (cuts.size() >= 30) return cuts;
                }
            }
        }
        return cuts;
    }
};

class BalasCutSeparator : public ICutSeparator {
public:
    std::vector<CutConstraint> separate(
        const std::vector<double> &primal_solution,
        const std::vector<double> &dual_solution,
        const BaseRelaxationModel &base,
        int ncols,
        double tol,
        double incumbent_bound) const override {
        std::vector<CutConstraint> cuts;
        if (static_cast<int>(dual_solution.size()) < base.nrows ||
            static_cast<int>(primal_solution.size()) < ncols)
            return cuts;

        std::vector<double> reduced_costs = compute_reduced_costs(
            base.obj, dual_solution, base, ncols);

        // Try Algorithm 2 (BCG) first
        CutConstraint bcg_cut = balas_cut_generate(
            primal_solution, dual_solution, reduced_costs,
            incumbent_bound, base, ncols, tol);

        if (!bcg_cut.indices.empty()) {
            // Check violation
            double lhs_val = 0.0;
            for (size_t k = 0; k < bcg_cut.indices.size(); ++k) {
                lhs_val += primal_solution[static_cast<size_t>(bcg_cut.indices[k])];
            }
            if (lhs_val < bcg_cut.rhs - tol) {
                cuts.push_back(std::move(bcg_cut));
            }
        }

        // Fall back to per-row heuristic for additional cuts
        struct CutWithViolation { CutConstraint cut; double violation; };
        std::vector<CutWithViolation> candidates;

        for (int i = 0; i < base.nrows; ++i) {
            if (dual_solution[static_cast<size_t>(i)] <= tol) continue;
            std::vector<int> Wi;
            double lhs_val = 0.0;
            for (int k = base.csr_offs[static_cast<size_t>(i)]; k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
                const int col = base.csr_inds[static_cast<size_t>(k)];
                if (col < 0 || col >= ncols) continue;
                if (base.csr_vals[static_cast<size_t>(k)] <= tol) continue;
                if (reduced_costs[static_cast<size_t>(col)] >= -tol) continue;
                if (primal_solution[static_cast<size_t>(col)] <= tol) continue;
                Wi.push_back(col);
                lhs_val += primal_solution[static_cast<size_t>(col)];
            }
            if (Wi.empty()) continue;
            const double violation = 1.0 - lhs_val;
            if (violation <= tol) continue;

            CutConstraint cut;
            cut.type = "balas_bcg";
            cut.rhs = 1.0;
            cut.indices = std::move(Wi);
            cut.values.assign(cut.indices.size(), 1.0);
            candidates.push_back({std::move(cut), violation});
        }

        std::sort(candidates.begin(), candidates.end(),
                  [](const CutWithViolation &a, const CutWithViolation &b) {
                      return a.violation > b.violation;
                  });

        const size_t max_cuts = 50;
        for (size_t ci = 0; ci < candidates.size() && cuts.size() < max_cuts; ++ci) {
            cuts.push_back(std::move(candidates[ci].cut));
        }
        return cuts;
    }
};

} // anonymous namespace

std::vector<std::unique_ptr<ICutSeparator>> make_cut_separators() {
    std::vector<std::unique_ptr<ICutSeparator>> separators;
    separators.push_back(std::make_unique<DualAggregatedCgSeparator>());
    separators.push_back(std::make_unique<RowPairCgSeparator>());
    separators.push_back(std::make_unique<BalasCutSeparator>());
    return separators;
}

void append_cuts_to_base_model(BaseRelaxationModel &base,
                                const std::vector<CutConstraint> &cuts) {
    if (cuts.empty()) return;

    for (const CutConstraint &cut : cuts) {
        for (size_t k = 0; k < cut.indices.size(); ++k) {
            base.csr_inds.push_back(cut.indices[k]);
            base.csr_vals.push_back(cut.values[k]);
        }
        base.csr_offs.push_back(static_cast<int>(base.csr_vals.size()));
        base.rhs.push_back(cut.rhs);
        base.base_cuts.push_back(cut);
    }

    base.nrows += static_cast<int>(cuts.size());
    base.nnz = static_cast<int>(base.csr_vals.size());
}

bool remap_cut_constraint(CutConstraint &cut, const std::vector<int> &old_to_new) {
    std::vector<int> new_indices;
    std::vector<double> new_values;
    new_indices.reserve(cut.indices.size());
    new_values.reserve(cut.values.size());

    for (size_t k = 0; k < cut.indices.size(); ++k) {
        const int old_col = cut.indices[k];
        if (old_col < 0 || old_col >= static_cast<int>(old_to_new.size())) continue;
        const int new_col = old_to_new[static_cast<size_t>(old_col)];
        if (new_col < 0) continue;
        new_indices.push_back(new_col);
        new_values.push_back(cut.values[k]);
    }

    cut.indices = std::move(new_indices);
    cut.values = std::move(new_values);
    return true;
}

} // namespace scpsol
