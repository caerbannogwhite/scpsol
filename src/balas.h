#ifndef SCPSOL_BALAS_H
#define SCPSOL_BALAS_H

#include "model.h"

namespace scpsol {

BalasBranchResult balas_branch_generate(
    const std::vector<double> &primal_solution,
    const std::vector<double> &dual_solution,
    const std::vector<double> &reduced_costs,
    const BaseRelaxationModel &base,
    int ncols,
    int max_branches,
    double integrality_tol,
    bool force_balas = false);

std::vector<BranchNodeState> balas_br1_children(
    const BranchNodeState &parent,
    const std::vector<std::vector<int>> &branch_sets,
    double parent_dual_bound,
    double parent_dual_bound_raw);

} // namespace scpsol

#endif // SCPSOL_BALAS_H
