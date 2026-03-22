#ifndef SCPSOL_MODEL_H
#define SCPSOL_MODEL_H

#include <chrono>
#include <limits>
#include <string>
#include <vector>

namespace scpsol {

struct ScpInstance {
    int nrows = 0;
    int ncols = 0;
    std::vector<double> costs;
    std::vector<int> csr_offsets;
    std::vector<int> csr_indices;
    std::vector<double> csr_values;
};

struct BranchDecision {
    int var_index;
    int fix_value; // 0 or 1
};

struct CutConstraint {
    std::vector<int> indices;
    std::vector<double> values;
    double rhs = 0.0;
    std::string type;
};

struct BranchNodeState {
    std::vector<BranchDecision> decisions;
    std::vector<CutConstraint> cuts;
    int depth = 0;
    double parent_dual_bound = -std::numeric_limits<double>::infinity();
    double parent_dual_bound_raw = -std::numeric_limits<double>::infinity();
};

struct BaseRelaxationModel {
    int nrows = 0;
    int ncols = 0;           // original variables only (no slacks)
    int ncols_input = 0;     // columns in input instance
    int nnz = 0;
    std::vector<int> csr_inds;
    std::vector<int> csr_offs;
    std::vector<double> csr_vals;
    std::vector<double> obj;
    std::vector<double> rhs;
    std::vector<int> active_to_original;
    std::vector<CutConstraint> base_cuts; // cuts appended to base model

    // Cached column-to-rows transpose
    struct ColRowEntry { int row; double val; };
    std::vector<std::vector<ColRowEntry>> cols_to_rows;

    void build_transpose() {
        cols_to_rows.assign(static_cast<size_t>(ncols), {});
        for (int i = 0; i < nrows; ++i) {
            for (int k = csr_offs[static_cast<size_t>(i)]; k < csr_offs[static_cast<size_t>(i) + 1]; ++k) {
                const int col = csr_inds[static_cast<size_t>(k)];
                if (col >= 0 && col < ncols)
                    cols_to_rows[static_cast<size_t>(col)].push_back({i, csr_vals[static_cast<size_t>(k)]});
            }
        }
    }
};

struct LpSolution {
    bool solved = false;
    bool optimal = false;
    bool infeasible = false;
    double primal_obj = std::numeric_limits<double>::infinity();
    double dual_obj = -std::numeric_limits<double>::infinity();
    std::vector<double> col_value;
    std::vector<double> row_dual;
    std::vector<double> col_dual;
};

struct IntegerHeuristicResult {
    bool feasible = false;
    double objective = std::numeric_limits<double>::infinity();
    std::vector<double> solution;
    std::string name;
};

struct BalasBranchResult {
    std::vector<std::vector<int>> sets;
    bool use_br1 = false;
};

struct ModelReductionResult {
    int columns_removed = 0;
    std::vector<int> old_to_new;
};

struct SolverConfig {
    int verbosity = 2;
    int max_nodes = 100000;
    double time_limit_seconds = 0.0;
    bool cuts_enabled = true;
    bool balas_enabled = true;
    int cut_rounds_root = 5;
    int max_cuts_per_round = 100;
    int balas_max_branches = 20;
    int gap_stagnation_window = 50;
    double mid_bnb_cut_frequency = 0.0;
    double aggressive_balas_frequency = 0.6;
    double lagrangian_frequency = 0.0; // 0.0 = disabled, 1.0 = every node
    double diving_frequency = 0.0; // 0.0 = disabled, 1.0 = every node
    double rins_frequency = 0.0;   // 0.0 = disabled, 1.0 = every node
    int diving_max_lp_solves = 50;
    std::string root_lp_solver = "auto"; // "simplex", "ipm", "auto"
    std::string decomposition_mode = "off"; // "off", "auto", "force"
    int decomposition_max_linkers = 20;
    int mid_bnb_cut_rounds = 3;
    int heuristic_frequency = 10;
    double integrality_tol = 1e-6;
    double feasibility_tol = 1e-6;
    double mip_gap_tol = 1e-6;
    double preprocess_time_limit = 10.0;
    bool show_solution = false;
    std::string preprocess_rules = "single,two";
    std::string branch_strategy = "reliability";
    int reliability_eta = 4;
    int reliability_max_sb = 8;
    std::string heuristic_config = "";
    double log_interval_seconds = 5.0;
};

struct SolverResult {
    std::string status;
    double primal_obj = std::numeric_limits<double>::infinity();
    double dual_obj = -std::numeric_limits<double>::infinity();
    double mip_gap = std::numeric_limits<double>::infinity();
    int nodes_processed = 0;
    int lp_solves = 0;
    double wall_time = 0.0;
    std::vector<double> solution;
};

struct ColumnPreprocessContext {
    int nrows = 0;
    int ncols = 0;
    std::vector<std::vector<int>> rows_by_column;
    std::vector<double> costs;
    std::vector<char> active;
    std::chrono::steady_clock::time_point deadline =
        std::chrono::steady_clock::time_point::max();
    double incumbent_bound = std::numeric_limits<double>::infinity();
};

} // namespace scpsol

#endif // SCPSOL_MODEL_H
