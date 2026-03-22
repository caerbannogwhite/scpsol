#ifndef SCPSOL_DIVING_H
#define SCPSOL_DIVING_H

#include "model.h"
#include "lp.h"

namespace scpsol {

struct DivingResult {
    bool found = false;
    double objective = std::numeric_limits<double>::infinity();
    std::vector<double> solution; // binary, ncols-dimensional
    int lp_solves = 0;
    std::string strategy_name;
};

// Run diving strategies at a node. Returns the best result found.
DivingResult run_diving_heuristics(
    LpSolver &lp,
    const BaseRelaxationModel &base,
    const BranchNodeState &branch_node,
    const LpSolution &current_sol,
    double best_obj,
    double integ_tol,
    int max_lp_solves);

} // namespace scpsol

#endif // SCPSOL_DIVING_H
