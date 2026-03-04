#ifndef SCPSOL_SOLVER_H
#define SCPSOL_SOLVER_H

#include "model.h"

namespace scpsol {

SolverResult solve(const ScpInstance &instance, const SolverConfig &config);

} // namespace scpsol

#endif // SCPSOL_SOLVER_H
