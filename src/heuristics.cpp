#include "heuristics.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <sstream>

namespace scpsol {

namespace {

class MostFractionalSelector : public IBranchVariableSelector {
public:
    int select(const std::vector<double> &solution,
               const std::vector<double> & /*objective*/,
               const std::vector<int> &candidates) const override {
        int best = -1;
        double best_score = -1.0;
        for (int idx : candidates) {
            const double frac = fabs(solution[static_cast<size_t>(idx)] -
                                     floor(solution[static_cast<size_t>(idx)] + 0.5));
            if (frac > best_score) {
                best_score = frac;
                best = idx;
            }
        }
        return best;
    }
};

class HighestCostFractionalSelector : public IBranchVariableSelector {
public:
    int select(const std::vector<double> & /*solution*/,
               const std::vector<double> &objective,
               const std::vector<int> &candidates) const override {
        int best = -1;
        double best_cost = -std::numeric_limits<double>::infinity();
        for (int idx : candidates) {
            if (objective[static_cast<size_t>(idx)] > best_cost) {
                best_cost = objective[static_cast<size_t>(idx)];
                best = idx;
            }
        }
        return best;
    }
};

class NearestIntegerFixingHeuristic : public IIntegerHeuristic {
public:
    IntegerHeuristicResult tryBuild(
        const std::vector<double> &relaxed_primal,
        const std::vector<double> & /*relaxed_dual*/,
        const BaseRelaxationModel &base,
        const BranchNodeState &branch_node,
        double tol) const override {
        IntegerHeuristicResult out;
        out.name = "nearest_integer_fixing";
        out.solution.assign(static_cast<size_t>(base.ncols), 0.0);

        if (static_cast<int>(relaxed_primal.size()) < base.ncols) return out;

        for (int j = 0; j < base.ncols; ++j) {
            const double rounded = floor(relaxed_primal[static_cast<size_t>(j)] + 0.5);
            out.solution[static_cast<size_t>(j)] = rounded < 0.0 ? 0.0 : (rounded > 1.0 ? 1.0 : rounded);
        }

        for (const BranchDecision &d : branch_node.decisions) {
            if (d.var_index >= 0 && d.var_index < base.ncols) {
                out.solution[static_cast<size_t>(d.var_index)] = static_cast<double>(d.fix_value);
            }
        }

        for (int i = 0; i < base.nrows; ++i) {
            double coverage = 0.0;
            for (int k = base.csr_offs[static_cast<size_t>(i)]; k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
                const int col = base.csr_inds[static_cast<size_t>(k)];
                if (col >= 0 && col < base.ncols) {
                    coverage += base.csr_vals[static_cast<size_t>(k)] * out.solution[static_cast<size_t>(col)];
                }
            }
            if (coverage + tol < base.rhs[static_cast<size_t>(i)]) return out;
        }

        out.feasible = true;
        out.objective = 0.0;
        for (int j = 0; j < base.ncols; ++j) {
            out.objective += base.obj[static_cast<size_t>(j)] * out.solution[static_cast<size_t>(j)];
        }
        return out;
    }
};

class DualGuidedCoverRepairHeuristic : public IIntegerHeuristic {
public:
    IntegerHeuristicResult tryBuild(
        const std::vector<double> &relaxed_primal,
        const std::vector<double> &relaxed_dual,
        const BaseRelaxationModel &base,
        const BranchNodeState &branch_node,
        double tol) const override {
        IntegerHeuristicResult out;
        out.name = "dual_guided_cover_repair";
        out.solution.assign(static_cast<size_t>(base.ncols), 0.0);

        if (static_cast<int>(relaxed_primal.size()) < base.ncols) return out;

        const int nrows = base.nrows;
        const int ncols = base.ncols;
        std::vector<char> fixed_zero(static_cast<size_t>(ncols), 0);
        std::vector<char> fixed_one(static_cast<size_t>(ncols), 0);
        std::vector<double> coverage(static_cast<size_t>(nrows), 0.0);

        // Use cached transpose from base model
        const auto &rows_by_col = base.cols_to_rows;
        if (static_cast<int>(rows_by_col.size()) < ncols) return out;

        for (const BranchDecision &d : branch_node.decisions) {
            if (d.var_index < 0 || d.var_index >= ncols) continue;
            if (d.fix_value == 0)
                fixed_zero[static_cast<size_t>(d.var_index)] = 1;
            else {
                fixed_one[static_cast<size_t>(d.var_index)] = 1;
                out.solution[static_cast<size_t>(d.var_index)] = 1.0;
            }
        }

        for (int j = 0; j < ncols; ++j) {
            if (fixed_zero[static_cast<size_t>(j)]) { out.solution[static_cast<size_t>(j)] = 0.0; continue; }
            if (fixed_one[static_cast<size_t>(j)]) continue;
            if (relaxed_primal[static_cast<size_t>(j)] >= 1.0 - tol) {
                out.solution[static_cast<size_t>(j)] = 1.0;
            }
        }

        // Compute initial coverage incrementally from selected columns
        for (int j = 0; j < ncols; ++j) {
            if (out.solution[static_cast<size_t>(j)] < 0.5) continue;
            for (const auto &entry : rows_by_col[static_cast<size_t>(j)])
                coverage[static_cast<size_t>(entry.row)] += entry.val;
        }

        auto is_row_covered = [&](int row) {
            return coverage[static_cast<size_t>(row)] + tol >= base.rhs[static_cast<size_t>(row)];
        };

        // Add coverage from column j
        auto add_col_coverage = [&](int j) {
            for (const auto &entry : rows_by_col[static_cast<size_t>(j)])
                coverage[static_cast<size_t>(entry.row)] += entry.val;
        };

        // Remove coverage from column j
        auto remove_col_coverage = [&](int j) {
            for (const auto &entry : rows_by_col[static_cast<size_t>(j)])
                coverage[static_cast<size_t>(entry.row)] -= entry.val;
        };

        while (true) {
            int uncovered = -1;
            for (int i = 0; i < nrows; ++i) {
                if (!is_row_covered(i)) { uncovered = i; break; }
            }
            if (uncovered < 0) break;

            // Score each free column using its transpose entries (O(nnz) total)
            int best_col = -1;
            double best_score = -std::numeric_limits<double>::infinity();
            for (int j = 0; j < ncols; ++j) {
                if (out.solution[static_cast<size_t>(j)] > 0.5 || fixed_zero[static_cast<size_t>(j)]) continue;
                double uncovered_gain = 0.0;
                double dual_gain = 0.0;
                for (const auto &entry : rows_by_col[static_cast<size_t>(j)]) {
                    if (is_row_covered(entry.row) || entry.val <= 0.0) continue;
                    uncovered_gain += entry.val;
                    if (entry.row < static_cast<int>(relaxed_dual.size()))
                        dual_gain += std::max(0.0, relaxed_dual[static_cast<size_t>(entry.row)]) * entry.val;
                }
                if (uncovered_gain <= 0.0) continue;
                const double col_cost = std::max(1e-9, base.obj[static_cast<size_t>(j)]);
                const double score = (uncovered_gain + dual_gain) / col_cost;
                if (score > best_score) { best_score = score; best_col = j; }
            }

            if (best_col < 0) {
                // Fallback: pick cheapest column covering any uncovered row
                int fallback_col = -1;
                double best_fallback_cost = std::numeric_limits<double>::infinity();
                for (int i = 0; i < nrows; ++i) {
                    if (is_row_covered(i)) continue;
                    for (int k = base.csr_offs[static_cast<size_t>(i)]; k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
                        const int col = base.csr_inds[static_cast<size_t>(k)];
                        if (col < 0 || col >= ncols || fixed_zero[static_cast<size_t>(col)] || out.solution[static_cast<size_t>(col)] > 0.5) continue;
                        if (base.csr_vals[static_cast<size_t>(k)] <= 0.0) continue;
                        if (base.obj[static_cast<size_t>(col)] < best_fallback_cost) {
                            best_fallback_cost = base.obj[static_cast<size_t>(col)];
                            fallback_col = col;
                        }
                    }
                }
                if (fallback_col < 0) return out;
                best_col = fallback_col;
            }
            out.solution[static_cast<size_t>(best_col)] = 1.0;
            add_col_coverage(best_col);
        }

        // Post-process: try removing unnecessary columns (most expensive first)
        std::vector<int> selected;
        for (int j = 0; j < ncols; ++j) {
            if (out.solution[static_cast<size_t>(j)] > 0.5 && !fixed_one[static_cast<size_t>(j)]) {
                selected.push_back(j);
            }
        }
        std::sort(selected.begin(), selected.end(),
                  [&](int a, int b) { return base.obj[static_cast<size_t>(a)] > base.obj[static_cast<size_t>(b)]; });

        for (int col : selected) {
            out.solution[static_cast<size_t>(col)] = 0.0;
            remove_col_coverage(col);
            bool feasible = true;
            for (int i = 0; i < nrows; ++i) {
                if (!is_row_covered(i)) { feasible = false; break; }
            }
            if (!feasible) {
                out.solution[static_cast<size_t>(col)] = 1.0;
                add_col_coverage(col);
            }
        }

        for (int i = 0; i < nrows; ++i) {
            if (!is_row_covered(i)) return out;
        }

        out.feasible = true;
        out.objective = 0.0;
        for (int j = 0; j < ncols; ++j) {
            out.objective += base.obj[static_cast<size_t>(j)] * out.solution[static_cast<size_t>(j)];
        }
        return out;
    }
};

std::string toLowerCopy(const std::string &s) {
    std::string out = s;
    std::transform(out.begin(), out.end(), out.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return out;
}

std::vector<std::string> splitCsvTokens(const std::string &csv) {
    std::vector<std::string> tokens;
    std::stringstream ss(csv);
    std::string item;
    while (std::getline(ss, item, ',')) {
        std::string cleaned;
        cleaned.reserve(item.size());
        for (char c : item) {
            if (!std::isspace(static_cast<unsigned char>(c))) cleaned.push_back(c);
        }
        if (!cleaned.empty()) tokens.push_back(toLowerCopy(cleaned));
    }
    return tokens;
}

} // anonymous namespace

std::unique_ptr<IBranchVariableSelector> make_branch_selector(const std::string &strategy) {
    const std::string s = toLowerCopy(strategy);
    if (s == "highest_cost_fractional")
        return std::make_unique<HighestCostFractionalSelector>();
    return std::make_unique<MostFractionalSelector>();
}

std::vector<std::unique_ptr<IIntegerHeuristic>> make_integer_heuristics(const std::string &configured) {
    std::vector<std::unique_ptr<IIntegerHeuristic>> heuristics;
    const auto tokens = splitCsvTokens(configured);
    if (tokens.empty()) {
        heuristics.push_back(std::make_unique<NearestIntegerFixingHeuristic>());
        heuristics.push_back(std::make_unique<DualGuidedCoverRepairHeuristic>());
        return heuristics;
    }
    for (const std::string &token : tokens) {
        if (token == "nearest_integer_fixing")
            heuristics.push_back(std::make_unique<NearestIntegerFixingHeuristic>());
        else if (token == "dual_guided_cover_repair")
            heuristics.push_back(std::make_unique<DualGuidedCoverRepairHeuristic>());
    }
    if (heuristics.empty()) {
        heuristics.push_back(std::make_unique<NearestIntegerFixingHeuristic>());
        heuristics.push_back(std::make_unique<DualGuidedCoverRepairHeuristic>());
    }
    return heuristics;
}

} // namespace scpsol
