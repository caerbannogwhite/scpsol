#ifndef SCPSOL_HEURISTICS_H
#define SCPSOL_HEURISTICS_H

#include "model.h"
#include <memory>

namespace scpsol {

class IBranchVariableSelector {
public:
    virtual ~IBranchVariableSelector() {}
    virtual int select(const std::vector<double> &solution,
                       const std::vector<double> &objective,
                       const std::vector<int> &candidates) const = 0;
};

class IIntegerHeuristic {
public:
    virtual ~IIntegerHeuristic() {}
    virtual IntegerHeuristicResult tryBuild(
        const std::vector<double> &relaxed_primal,
        const std::vector<double> &relaxed_dual,
        const BaseRelaxationModel &base,
        const BranchNodeState &branch_node,
        double tol) const = 0;
};

std::unique_ptr<IBranchVariableSelector> make_branch_selector(const std::string &strategy);
std::vector<std::unique_ptr<IIntegerHeuristic>> make_integer_heuristics(const std::string &configured);

} // namespace scpsol

#endif // SCPSOL_HEURISTICS_H
