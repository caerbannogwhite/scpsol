#ifndef SCPSOL_RELIABILITY_H
#define SCPSOL_RELIABILITY_H

#include "model.h"
#include "lp.h"

namespace scpsol {

struct PseudocostState {
    std::vector<double> sum_down;   // sum of (gain / frac) per variable
    std::vector<double> sum_up;     // sum of (gain / (1-frac)) per variable
    std::vector<int> count_down;
    std::vector<int> count_up;

    void init(int ncols);
    void update_down(int j, double frac, double gain);
    void update_up(int j, double frac, double gain);
    double estimate_down(int j, double frac, double obj_j) const;
    double estimate_up(int j, double frac, double obj_j) const;
    bool is_reliable(int j, int eta_rel) const;
    void remap(const std::vector<int> &old_to_new, int new_ncols);
};

// Lightweight propagation: essential-column cascade from a set of fixings.
struct PropagationResult {
    std::vector<BranchDecision> implied;
    bool infeasible = false;
};

// Pre-built adjacency for propagation (built once per node).
struct ScpAdjacency {
    std::vector<std::vector<int>> cols_by_row;
    std::vector<std::vector<int>> rows_by_col;
    int nrows = 0;
    int ncols = 0;
};

ScpAdjacency build_adjacency(const BaseRelaxationModel &base);

PropagationResult propagate_fixings(
    const ScpAdjacency &adj,
    const std::vector<BranchDecision> &current_decisions,
    int probe_var, int probe_value);

// Reliability branching variable selection.
// Returns variable index to branch on, or -1 if no candidate.
// Updates pc_state when strong branching is performed.
// Increments lp_solves_out by the number of probe LPs solved.
int reliability_branch_select(
    PseudocostState &pc_state,
    LpSolver &lp,
    const BaseRelaxationModel &base,
    const BranchNodeState &branch_node,
    const LpSolution &current_sol,
    const std::vector<int> &fractional_candidates,
    double best_obj,
    int eta_rel,
    int max_sb_candidates,
    double integ_tol,
    int verbosity,
    int &lp_solves_out,
    const ScpAdjacency *cached_adj = nullptr);

} // namespace scpsol

#endif // SCPSOL_RELIABILITY_H
