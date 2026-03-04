#ifndef SCPSOL_PREPROCESSOR_H
#define SCPSOL_PREPROCESSOR_H

#include "model.h"
#include <memory>

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

class IColumnPreprocessRule {
public:
    virtual ~IColumnPreprocessRule() {}
    virtual const char *name() const = 0;
    virtual int apply(ColumnPreprocessContext &ctx, double tol) const = 0;
};

std::vector<std::unique_ptr<IColumnPreprocessRule>>
make_preprocess_rules(const std::string &configured);

} // namespace scpsol

#endif // SCPSOL_PREPROCESSOR_H
