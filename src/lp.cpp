#include "lp.h"

#include <algorithm>
#include <cstdio>
#include <limits>

namespace scpsol {

static void csr_to_csc(int nrows, int ncols,
                        const std::vector<int> &csr_offs,
                        const std::vector<int> &csr_inds,
                        const std::vector<double> &csr_vals,
                        std::vector<int> &csc_start,
                        std::vector<int> &csc_index,
                        std::vector<double> &csc_value) {
    const int nnz = static_cast<int>(csr_vals.size());
    csc_start.assign(static_cast<size_t>(ncols + 1), 0);
    csc_index.resize(static_cast<size_t>(nnz));
    csc_value.resize(static_cast<size_t>(nnz));

    // Count entries per column
    for (int k = 0; k < nnz; ++k) {
        ++csc_start[static_cast<size_t>(csr_inds[static_cast<size_t>(k)] + 1)];
    }
    // Prefix sum
    for (int j = 1; j <= ncols; ++j) {
        csc_start[static_cast<size_t>(j)] += csc_start[static_cast<size_t>(j - 1)];
    }
    // Fill
    std::vector<int> cursor(csc_start.begin(), csc_start.begin() + ncols);
    for (int i = 0; i < nrows; ++i) {
        for (int k = csr_offs[static_cast<size_t>(i)]; k < csr_offs[static_cast<size_t>(i + 1)]; ++k) {
            const int col = csr_inds[static_cast<size_t>(k)];
            const int pos = cursor[static_cast<size_t>(col)]++;
            csc_index[static_cast<size_t>(pos)] = i;
            csc_value[static_cast<size_t>(pos)] = csr_vals[static_cast<size_t>(k)];
        }
    }
}

void LpSolver::build_model(const BaseRelaxationModel &base) {
    highs_.clear();
    highs_.setOptionValue("output_flag", false);
    highs_.setOptionValue("presolve", "off");
    highs_.setOptionValue("simplex_strategy", 1); // dual simplex

    const int ncols = base.ncols;
    const int nrows = base.nrows;

    // Convert CSR to CSC for HiGHS
    std::vector<int> csc_start, csc_index;
    std::vector<double> csc_value;
    csr_to_csc(nrows, ncols, base.csr_offs, base.csr_inds, base.csr_vals,
               csc_start, csc_index, csc_value);

    // Column bounds: 0 <= x_j <= 1 for all original variables
    std::vector<double> col_lower(static_cast<size_t>(ncols), 0.0);
    std::vector<double> col_upper(static_cast<size_t>(ncols), 1.0);

    // Row bounds: >= 1.0 for cover constraints, >= rhs for cuts
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> row_lower(static_cast<size_t>(nrows));
    std::vector<double> row_upper(static_cast<size_t>(nrows), inf);
    for (int i = 0; i < nrows; ++i) {
        row_lower[static_cast<size_t>(i)] = base.rhs[static_cast<size_t>(i)];
    }

    HighsModel model;
    model.lp_.num_col_ = ncols;
    model.lp_.num_row_ = nrows;
    model.lp_.sense_ = ObjSense::kMinimize;
    model.lp_.col_cost_ = base.obj;
    model.lp_.col_lower_ = col_lower;
    model.lp_.col_upper_ = col_upper;
    model.lp_.row_lower_ = row_lower;
    model.lp_.row_upper_ = row_upper;
    model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
    model.lp_.a_matrix_.start_ = csc_start;
    model.lp_.a_matrix_.index_ = csc_index;
    model.lp_.a_matrix_.value_ = csc_value;

    highs_.passModel(model);

    base_num_rows_ = nrows;
    base_num_cols_ = ncols;
    base_lower_ = col_lower;
    base_upper_ = col_upper;
    has_saved_basis_ = false;
}

void LpSolver::apply_decisions(const std::vector<BranchDecision> &decisions) {
    for (const BranchDecision &d : decisions) {
        if (d.var_index < 0 || d.var_index >= base_num_cols_) continue;
        double lb = d.fix_value == 1 ? 1.0 : 0.0;
        double ub = d.fix_value == 0 ? 0.0 : 1.0;
        highs_.changeColBounds(d.var_index, lb, ub);
    }
}

void LpSolver::add_cuts(const std::vector<CutConstraint> &cuts) {
    const double inf = std::numeric_limits<double>::infinity();
    for (const CutConstraint &cut : cuts) {
        if (cut.indices.empty()) continue;
        highs_.addRow(cut.rhs, inf,
                      static_cast<HighsInt>(cut.indices.size()),
                      cut.indices.data(), cut.values.data());
    }
}

void LpSolver::restore_base_state() {
    // Remove any rows added beyond the base model
    const int current_rows = highs_.getNumRow();
    if (current_rows > base_num_rows_) {
        std::vector<HighsInt> rows_to_delete;
        for (int i = base_num_rows_; i < current_rows; ++i) {
            rows_to_delete.push_back(i);
        }
        highs_.deleteRows(static_cast<HighsInt>(rows_to_delete.size()),
                          rows_to_delete.data());
    }
    // Reset column bounds
    for (int j = 0; j < base_num_cols_; ++j) {
        highs_.changeColBounds(j, base_lower_[static_cast<size_t>(j)],
                               base_upper_[static_cast<size_t>(j)]);
    }
}

LpSolution LpSolver::solve() {
    LpSolution sol;

    if (has_saved_basis_) {
        highs_.setBasis(saved_basis_);
    }

    HighsStatus status = highs_.run();
    const HighsInfo &info = highs_.getInfo();
    const HighsModelStatus model_status = highs_.getModelStatus();

    sol.solved = (status == HighsStatus::kOk);
    sol.optimal = (model_status == HighsModelStatus::kOptimal);
    sol.infeasible = (model_status == HighsModelStatus::kInfeasible);

    if (sol.solved && sol.optimal) {
        sol.primal_obj = info.objective_function_value;
        sol.dual_obj = info.objective_function_value; // At optimality they match

        const HighsSolution &hs = highs_.getSolution();
        sol.col_value = hs.col_value;
        sol.row_dual = hs.row_dual;
        sol.col_dual = hs.col_dual;
    }

    return sol;
}

void LpSolver::save_basis() {
    saved_basis_ = highs_.getBasis();
    has_saved_basis_ = true;
}

void LpSolver::restore_basis() {
    if (has_saved_basis_) {
        highs_.setBasis(saved_basis_);
    }
}

void LpSolver::rebuild_model(const BaseRelaxationModel &base) {
    build_model(base);
}

void LpSolver::rebuild_model_keep_basis(const BaseRelaxationModel &base) {
    HighsBasis basis;
    bool had_basis = has_saved_basis_;
    if (had_basis) basis = saved_basis_;
    build_model(base);
    if (had_basis) {
        // Truncate or extend basis to match the new model dimensions
        basis.col_status.resize(static_cast<size_t>(base.ncols), HighsBasisStatus::kNonbasic);
        basis.row_status.resize(static_cast<size_t>(base.nrows), HighsBasisStatus::kNonbasic);
        saved_basis_ = basis;
        has_saved_basis_ = true;
    }
}

} // namespace scpsol
