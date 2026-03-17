#ifndef SCPSOL_LP_H
#define SCPSOL_LP_H

#include "model.h"
#include <Highs.h>

namespace scpsol {

class LpSolver {
public:
    void build_model(const BaseRelaxationModel &base);
    void apply_decisions(const std::vector<BranchDecision> &decisions);
    void add_cuts(const std::vector<CutConstraint> &cuts);
    void restore_base_state();
    LpSolution solve();
    void save_basis();
    void restore_basis();
    HighsBasis get_basis() const;
    void set_basis(const HighsBasis &basis);
    void rebuild_model(const BaseRelaxationModel &base);
    void rebuild_model_keep_basis(const BaseRelaxationModel &base);

private:
    Highs highs_;
    int base_num_rows_ = 0;
    int base_num_cols_ = 0;
    HighsBasis saved_basis_;
    bool has_saved_basis_ = false;
    std::vector<double> base_lower_;
    std::vector<double> base_upper_;
};

} // namespace scpsol

#endif // SCPSOL_LP_H
