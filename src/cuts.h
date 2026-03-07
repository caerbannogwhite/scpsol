#ifndef SCPSOL_CUTS_H
#define SCPSOL_CUTS_H

#include "model.h"
#include <memory>

namespace scpsol {

class ICutSeparator {
public:
    virtual ~ICutSeparator() {}
    virtual std::vector<CutConstraint> separate(
        const std::vector<double> &primal_solution,
        const std::vector<double> &dual_solution,
        const BaseRelaxationModel &base,
        int ncols,
        double tol,
        double incumbent_bound) const = 0;
};

std::vector<std::unique_ptr<ICutSeparator>> make_cut_separators();

void append_cuts_to_base_model(BaseRelaxationModel &base,
                                const std::vector<CutConstraint> &cuts);

bool remap_cut_constraint(CutConstraint &cut, const std::vector<int> &old_to_new);

std::vector<double> compute_reduced_costs(
    const std::vector<double> &objective,
    const std::vector<double> &dual_solution,
    const BaseRelaxationModel &base,
    int ncols);

} // namespace scpsol

#endif // SCPSOL_CUTS_H
