#ifndef SCPSOL_DECOMPOSITION_H
#define SCPSOL_DECOMPOSITION_H

#include "model.h"
#include <vector>

namespace scpsol {

struct ScpAdjacency; // forward (defined in reliability.h)

struct DecompositionState {
    bool enabled = false;
    bool computed = false;
    std::vector<int> linking_cols; // indices in base model space

    struct BlockInfo {
        std::vector<int> rows; // row indices in parent model
        std::vector<int> cols; // col indices in parent model
    };
    std::vector<BlockInfo> blocks;
};

// Find a small set of linking columns whose removal disconnects the
// residual constraint matrix into independent blocks.
DecompositionState compute_decomposition(
    const BaseRelaxationModel &base,
    int max_linkers,
    int verbosity);

// Check whether all linking columns are fixed in the given decisions.
bool all_linkers_fixed(
    const DecompositionState &state,
    const std::vector<BranchDecision> &decisions);

// Pick the best unfixed linking column for priority branching.
// Returns -1 if no unfixed linking column is among the fractional candidates.
int pick_linking_branch_var(
    const DecompositionState &state,
    const std::vector<BranchDecision> &decisions,
    const std::vector<int> &fractional);

struct DecompositionResult {
    bool solved = false;        // true if all sub-problems were solved
    bool improved = false;      // true if a better incumbent was found
    double combined_obj = std::numeric_limits<double>::infinity();
    std::vector<double> combined_solution; // in parent model column space
};

// Solve the decomposed sub-problems once all linking columns are fixed.
DecompositionResult solve_decomposed(
    const BaseRelaxationModel &base,
    const std::vector<BranchDecision> &decisions,
    const DecompositionState &state,
    double incumbent_obj,
    const SolverConfig &config,
    double time_remaining,
    int verbosity);

} // namespace scpsol

#endif // SCPSOL_DECOMPOSITION_H
