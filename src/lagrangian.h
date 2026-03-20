#ifndef SCPSOL_LAGRANGIAN_H
#define SCPSOL_LAGRANGIAN_H

#include "model.h"
#include <chrono>
#include <vector>

namespace scpsol {

struct LagrangianResult {
    double best_bound = -std::numeric_limits<double>::infinity();
    double best_heuristic_obj = std::numeric_limits<double>::infinity();
    std::vector<double> best_solution; // binary feasible solution (ncols)
    std::vector<double> best_multipliers;
    int iterations = 0;
};

LagrangianResult lagrangian_relaxation(
    const BaseRelaxationModel &base,
    double incumbent_obj,
    const std::vector<double> &init_multipliers,
    int max_iterations,
    double time_limit,
    int verbosity);

// Node-level variant: respects branch fixings from decisions.
LagrangianResult lagrangian_relaxation_at_node(
    const BaseRelaxationModel &base,
    double incumbent_obj,
    const std::vector<double> &init_multipliers,
    const std::vector<BranchDecision> &decisions,
    int max_iterations,
    double time_limit);

} // namespace scpsol

#endif // SCPSOL_LAGRANGIAN_H
