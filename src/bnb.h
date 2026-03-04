#ifndef SCPSOL_BNB_H
#define SCPSOL_BNB_H

#include "model.h"

namespace scpsol {

ModelReductionResult reduce_base_model(
    BaseRelaxationModel &base, double incumbent_bound, double tol);

ModelReductionResult reduce_base_model_budget_pruning(
    BaseRelaxationModel &base, double incumbent_bound, double tol,
    double preprocess_time_limit_sec);

bool remap_branch_node(BranchNodeState &node, const std::vector<int> &old_to_new);

bool append_decision_if_consistent(
    const BranchNodeState &parent, int var, int value, BranchNodeState *child);

bool is_binary_integral(const std::vector<double> &x, int ncols, double tol);

std::vector<int> collect_fractional_candidates(
    const std::vector<double> &x, int ncols, double tol);

bool is_node_provably_infeasible(
    const BranchNodeState &node, const BaseRelaxationModel &base);

bool has_integer_objective(const std::vector<double> &obj, int ncols, double tol);

double tighten_dual_bound(double bound, double tol);

double compute_mip_gap(double incumbent, double dual_bound);

} // namespace scpsol

#endif // SCPSOL_BNB_H
