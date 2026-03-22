#ifndef SCPSOL_RINS_H
#define SCPSOL_RINS_H

#include "model.h"

namespace scpsol {

struct RinsResult {
    bool found = false;
    double objective = std::numeric_limits<double>::infinity();
    std::vector<double> solution; // in parent model's active column space
    int sub_nodes = 0;
    int sub_lp_solves = 0;
};

// RINS: fix variables where LP relaxation and incumbent agree,
// solve the residual SCP on the disagreeing variables.
RinsResult run_rins(
    const BaseRelaxationModel &base,
    const LpSolution &lp_sol,
    const std::vector<double> &incumbent_active, // incumbent in active col space
    double incumbent_obj,
    const SolverConfig &config,
    double time_budget,
    int verbosity);

} // namespace scpsol

#endif // SCPSOL_RINS_H
