#ifndef SCPSOL_PREPROCESSOR_H
#define SCPSOL_PREPROCESSOR_H

#include "model.h"
#include <memory>
#include <string>

namespace scpsol {

struct GreedySetCoverResult {
    bool feasible = false;
    double objective = 0.0;
    std::vector<int> selected_columns;
};

GreedySetCoverResult greedy_set_cover_heuristic(
    int nrows, int ncols,
    const std::vector<int> &csr_inds,
    const std::vector<int> &csr_offs,
    const std::vector<double> &csr_vals,
    const double *obj);

struct RowReductionResult {
    int rows_removed = 0;
    int cols_fixed = 0;
    double fixed_cost = 0.0;
    std::vector<int> fixed_original_cols;
};

RowReductionResult row_reduce(
    int &nrows, int &ncols,
    std::vector<int> &csr_inds, std::vector<int> &csr_offs, std::vector<double> &csr_vals,
    std::vector<double> &obj, std::vector<int> &active_to_input,
    double tol, double time_limit_sec, int verbosity);

// Build the column-oriented view (rows per column, costs, all-active flags,
// optional deadline and incumbent bound) that the column rules operate on.
// A time limit <= 0 means no deadline.
ColumnPreprocessContext make_column_context(
    int nrows, int ncols,
    const std::vector<int> &csr_inds, const std::vector<int> &csr_offs,
    const std::vector<double> &csr_vals, const std::vector<double> &obj,
    double tol, double time_limit_sec, double incumbent_bound);

class IColumnPreprocessRule {
public:
    virtual ~IColumnPreprocessRule() {}
    virtual const char *name() const = 0;
    virtual int apply(ColumnPreprocessContext &ctx, double tol) const = 0;
};

std::vector<std::unique_ptr<IColumnPreprocessRule>>
make_preprocess_rules(const std::string &configured);

// Dominance Finder: for each non-unit-cost column, use greedy set cover to
// check if cheaper columns can collectively cover the same rows at lower cost.
// Returns number of columns removed. Modifies ctx.active in place.
int dominance_finder(ColumnPreprocessContext &ctx, double tol,
                     double time_limit_per_sub_mip, int verbosity);

} // namespace scpsol

#endif // SCPSOL_PREPROCESSOR_H
