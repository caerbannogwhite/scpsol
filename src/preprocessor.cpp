#include "preprocessor.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <sstream>

namespace scpsol {

GreedySetCoverResult greedy_set_cover_heuristic(
    int nrows, int ncols,
    const std::vector<int> &csr_inds,
    const std::vector<int> &csr_offs,
    const std::vector<double> &csr_vals,
    const double *obj) {
    GreedySetCoverResult result;
    if (nrows <= 0 || ncols <= 0) return result;

    // Build rowsByColumn (transpose CSR)
    std::vector<std::vector<int>> rows_by_col(static_cast<size_t>(ncols));
    for (int i = 0; i < nrows; ++i) {
        const int begin = csr_offs[static_cast<size_t>(i)];
        const int end = csr_offs[static_cast<size_t>(i) + 1];
        for (int k = begin; k < end; ++k) {
            const int col = csr_inds[static_cast<size_t>(k)];
            if (col >= 0 && col < ncols && csr_vals[static_cast<size_t>(k)] > 0.0) {
                rows_by_col[static_cast<size_t>(col)].push_back(i);
            }
        }
    }

    struct ColEntry { double cost; int neg_coverage; int col_index; };
    std::vector<ColEntry> columns(static_cast<size_t>(ncols));
    for (int j = 0; j < ncols; ++j) {
        columns[static_cast<size_t>(j)] = {
            obj[j],
            -static_cast<int>(rows_by_col[static_cast<size_t>(j)].size()),
            j
        };
    }
    std::sort(columns.begin(), columns.end(), [](const ColEntry &a, const ColEntry &b) {
        if (a.cost != b.cost) return a.cost < b.cost;
        return a.neg_coverage < b.neg_coverage;
    });

    std::vector<char> covered(static_cast<size_t>(nrows), 0);
    int uncovered_count = nrows;
    double total_cost = 0.0;

    for (const ColEntry &entry : columns) {
        if (uncovered_count <= 0) break;
        int new_coverage = 0;
        for (int row : rows_by_col[static_cast<size_t>(entry.col_index)]) {
            if (!covered[static_cast<size_t>(row)]) ++new_coverage;
        }
        if (new_coverage > 0) {
            for (int row : rows_by_col[static_cast<size_t>(entry.col_index)]) {
                if (!covered[static_cast<size_t>(row)]) {
                    covered[static_cast<size_t>(row)] = 1;
                    --uncovered_count;
                }
            }
            total_cost += entry.cost;
            result.selected_columns.push_back(entry.col_index);
        }
    }

    if (uncovered_count == 0) {
        result.feasible = true;
        result.objective = total_cost;
    }
    return result;
}

namespace {

bool isSubsetSorted(const std::vector<int> &subset, const std::vector<int> &superset) {
    size_t i = 0, j = 0;
    while (i < subset.size() && j < superset.size()) {
        if (subset[i] == superset[j]) { ++i; ++j; }
        else if (subset[i] > superset[j]) { ++j; }
        else { return false; }
    }
    return i == subset.size();
}

bool unionCoversSorted(const std::vector<int> &target,
                       const std::vector<int> &a,
                       const std::vector<int> &b) {
    size_t ia = 0, ib = 0;
    for (int row : target) {
        while (ia < a.size() && a[ia] < row) ++ia;
        if (ia < a.size() && a[ia] == row) continue;
        while (ib < b.size() && b[ib] < row) ++ib;
        if (ib < b.size() && b[ib] == row) continue;
        return false;
    }
    return true;
}

bool tripleUnionCoversSorted(const std::vector<int> &target,
                              const std::vector<int> &a,
                              const std::vector<int> &b,
                              const std::vector<int> &c) {
    size_t ia = 0, ib = 0, ic = 0;
    for (int row : target) {
        while (ia < a.size() && a[ia] < row) ++ia;
        if (ia < a.size() && a[ia] == row) continue;
        while (ib < b.size() && b[ib] < row) ++ib;
        if (ib < b.size() && b[ib] == row) continue;
        while (ic < c.size() && c[ic] < row) ++ic;
        if (ic < c.size() && c[ic] == row) continue;
        return false;
    }
    return true;
}

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

class SingleColumnDominanceRule : public IColumnPreprocessRule {
public:
    const char *name() const override { return "single_column_dominance"; }
    int apply(ColumnPreprocessContext &ctx, double tol) const override {
        int removed = 0;
        for (int target = 0; target < ctx.ncols; ++target) {
            if (std::chrono::steady_clock::now() >= ctx.deadline) break;
            if (!ctx.active[static_cast<size_t>(target)]) continue;
            const auto &target_rows = ctx.rows_by_column[static_cast<size_t>(target)];
            for (int cand = 0; cand < ctx.ncols; ++cand) {
                if (cand == target || !ctx.active[static_cast<size_t>(cand)]) continue;
                if (ctx.costs[static_cast<size_t>(cand)] > ctx.costs[static_cast<size_t>(target)] + tol) continue;
                const auto &cand_rows = ctx.rows_by_column[static_cast<size_t>(cand)];
                if (!isSubsetSorted(target_rows, cand_rows)) continue;
                if (fabs(ctx.costs[static_cast<size_t>(cand)] - ctx.costs[static_cast<size_t>(target)]) <= tol && cand > target) continue;
                ctx.active[static_cast<size_t>(target)] = 0;
                ++removed;
                break;
            }
        }
        return removed;
    }
};

class TwoColumnDominanceRule : public IColumnPreprocessRule {
public:
    const char *name() const override { return "two_column_dominance"; }
    int apply(ColumnPreprocessContext &ctx, double tol) const override {
        int removed = 0;
        bool timed_out = false;
        for (int target = 0; target < ctx.ncols && !timed_out; ++target) {
            if (std::chrono::steady_clock::now() >= ctx.deadline) break;
            if (!ctx.active[static_cast<size_t>(target)]) continue;
            const auto &target_rows = ctx.rows_by_column[static_cast<size_t>(target)];
            const double target_cost = ctx.costs[static_cast<size_t>(target)];
            bool dominated = false;
            for (int a = 0; a < ctx.ncols && !dominated; ++a) {
                if ((a & 255) == 0 && std::chrono::steady_clock::now() >= ctx.deadline) { timed_out = true; break; }
                if (a == target || !ctx.active[static_cast<size_t>(a)]) continue;
                const double a_cost = ctx.costs[static_cast<size_t>(a)];
                if (a_cost >= target_cost - tol) continue;
                for (int b = a + 1; b < ctx.ncols; ++b) {
                    if (b == target || !ctx.active[static_cast<size_t>(b)]) continue;
                    const double pair_cost = a_cost + ctx.costs[static_cast<size_t>(b)];
                    if (pair_cost >= target_cost - tol) continue;
                    if (unionCoversSorted(target_rows, ctx.rows_by_column[static_cast<size_t>(a)],
                                          ctx.rows_by_column[static_cast<size_t>(b)])) {
                        dominated = true;
                        break;
                    }
                }
            }
            if (dominated) {
                ctx.active[static_cast<size_t>(target)] = 0;
                ++removed;
            }
        }
        return removed;
    }
};

class CostDrivenReplacementRule : public IColumnPreprocessRule {
public:
    const char *name() const override { return "cost_driven_replacement"; }
    int apply(ColumnPreprocessContext &ctx, double tol) const override {
        if (ctx.nrows <= 0 || ctx.ncols <= 0) return 0;

        std::vector<std::vector<int>> columns_by_row(static_cast<size_t>(ctx.nrows));
        for (int j = 0; j < ctx.ncols; ++j) {
            if (!ctx.active[static_cast<size_t>(j)]) continue;
            for (int row : ctx.rows_by_column[static_cast<size_t>(j)])
                columns_by_row[static_cast<size_t>(row)].push_back(j);
        }

        std::vector<int> sorted_cols;
        sorted_cols.reserve(static_cast<size_t>(ctx.ncols));
        for (int j = 0; j < ctx.ncols; ++j) {
            if (ctx.active[static_cast<size_t>(j)]) sorted_cols.push_back(j);
        }
        std::sort(sorted_cols.begin(), sorted_cols.end(),
                  [&](int a, int b) { return ctx.costs[static_cast<size_t>(a)] > ctx.costs[static_cast<size_t>(b)]; });

        int removed = 0;
        std::vector<char> seen(static_cast<size_t>(ctx.ncols), 0);
        bool timed_out = false;

        for (int target : sorted_cols) {
            if (std::chrono::steady_clock::now() >= ctx.deadline) break;
            if (!ctx.active[static_cast<size_t>(target)]) continue;
            const auto &target_rows = ctx.rows_by_column[static_cast<size_t>(target)];
            if (target_rows.empty()) continue;
            const double target_cost = ctx.costs[static_cast<size_t>(target)];

            std::vector<int> candidates;
            for (int row : target_rows) {
                for (int col : columns_by_row[static_cast<size_t>(row)]) {
                    if (col != target && ctx.active[static_cast<size_t>(col)] && !seen[static_cast<size_t>(col)]) {
                        seen[static_cast<size_t>(col)] = 1;
                        candidates.push_back(col);
                    }
                }
            }
            for (int col : candidates) seen[static_cast<size_t>(col)] = 0;

            std::sort(candidates.begin(), candidates.end(),
                      [&](int a, int b) { return ctx.costs[static_cast<size_t>(a)] < ctx.costs[static_cast<size_t>(b)]; });

            bool dominated = false;

            // Try pairs
            for (size_t i = 0; i < candidates.size() && !dominated; ++i) {
                if ((i & 255) == 0 && std::chrono::steady_clock::now() >= ctx.deadline) { timed_out = true; break; }
                const int a = candidates[i];
                const double cost_a = ctx.costs[static_cast<size_t>(a)];
                if (cost_a > target_cost + tol) break;
                for (size_t j = i + 1; j < candidates.size() && !dominated; ++j) {
                    const int b = candidates[j];
                    if (cost_a + ctx.costs[static_cast<size_t>(b)] > target_cost + tol) break;
                    if (unionCoversSorted(target_rows, ctx.rows_by_column[static_cast<size_t>(a)],
                                          ctx.rows_by_column[static_cast<size_t>(b)])) {
                        dominated = true;
                    }
                }
            }
            if (timed_out) break;

            // Try triplets
            if (!dominated) {
                for (size_t i = 0; i < candidates.size() && !dominated; ++i) {
                    if ((i & 63) == 0 && std::chrono::steady_clock::now() >= ctx.deadline) { timed_out = true; break; }
                    const int a = candidates[i];
                    const double cost_a = ctx.costs[static_cast<size_t>(a)];
                    if (cost_a > target_cost + tol) break;
                    for (size_t j = i + 1; j < candidates.size() && !dominated; ++j) {
                        const int b = candidates[j];
                        const double cost_ab = cost_a + ctx.costs[static_cast<size_t>(b)];
                        if (cost_ab > target_cost + tol) break;
                        for (size_t k = j + 1; k < candidates.size() && !dominated; ++k) {
                            const int c = candidates[k];
                            if (cost_ab + ctx.costs[static_cast<size_t>(c)] > target_cost + tol) break;
                            if (tripleUnionCoversSorted(target_rows, ctx.rows_by_column[static_cast<size_t>(a)],
                                                        ctx.rows_by_column[static_cast<size_t>(b)],
                                                        ctx.rows_by_column[static_cast<size_t>(c)])) {
                                dominated = true;
                            }
                        }
                    }
                }
            }
            if (timed_out) break;
            if (dominated) {
                ctx.active[static_cast<size_t>(target)] = 0;
                ++removed;
            }
        }
        return removed;
    }
};

class IncumbentBudgetPruningRule : public IColumnPreprocessRule {
public:
    const char *name() const override { return "incumbent_budget_pruning"; }
    int apply(ColumnPreprocessContext &ctx, double tol) const override {
        if (ctx.nrows <= 0 || ctx.ncols <= 0 || !std::isfinite(ctx.incumbent_bound)) return 0;

        const double incumbent_floor = std::floor(ctx.incumbent_bound);

        struct CostCol { double cost; int col; };
        std::vector<std::vector<CostCol>> columns_by_row(static_cast<size_t>(ctx.nrows));
        for (int j = 0; j < ctx.ncols; ++j) {
            if (!ctx.active[static_cast<size_t>(j)]) continue;
            for (int row : ctx.rows_by_column[static_cast<size_t>(j)]) {
                columns_by_row[static_cast<size_t>(row)].push_back({ctx.costs[static_cast<size_t>(j)], j});
            }
        }
        for (auto &vec : columns_by_row) {
            std::sort(vec.begin(), vec.end(), [](const CostCol &a, const CostCol &b) { return a.cost < b.cost; });
        }

        std::vector<int> sorted_cols;
        sorted_cols.reserve(static_cast<size_t>(ctx.ncols));
        for (int j = 0; j < ctx.ncols; ++j) {
            if (ctx.active[static_cast<size_t>(j)]) sorted_cols.push_back(j);
        }
        std::sort(sorted_cols.begin(), sorted_cols.end(),
                  [&](int a, int b) { return ctx.costs[static_cast<size_t>(a)] > ctx.costs[static_cast<size_t>(b)]; });

        std::vector<int> cost1_cols;
        for (int j = 0; j < ctx.ncols; ++j) {
            if (ctx.active[static_cast<size_t>(j)] &&
                std::fabs(ctx.costs[static_cast<size_t>(j)] - 1.0) <= tol) {
                cost1_cols.push_back(j);
            }
        }

        int removed = 0;
        for (int target : sorted_cols) {
            if (std::chrono::steady_clock::now() >= ctx.deadline) break;
            if (!ctx.active[static_cast<size_t>(target)]) continue;
            const double cost_j = ctx.costs[static_cast<size_t>(target)];
            const double budget = incumbent_floor - 1.0 - std::floor(cost_j);

            if (budget < -tol) {
                ctx.active[static_cast<size_t>(target)] = 0;
                ++removed;
                continue;
            }

            const auto &target_rows = ctx.rows_by_column[static_cast<size_t>(target)];
            std::vector<int> uncovered_rows;
            {
                size_t ti = 0;
                for (int r = 0; r < ctx.nrows; ++r) {
                    if (ti < target_rows.size() && target_rows[ti] == r) { ++ti; continue; }
                    uncovered_rows.push_back(r);
                }
            }

            if (uncovered_rows.empty()) continue;

            if (budget < tol) {
                ctx.active[static_cast<size_t>(target)] = 0;
                ++removed;
                continue;
            }

            if (budget < 1.0 + tol) {
                bool found = false;
                for (int k : cost1_cols) {
                    if (std::chrono::steady_clock::now() >= ctx.deadline) goto done;
                    if (k == target || !ctx.active[static_cast<size_t>(k)]) continue;
                    if (isSubsetSorted(uncovered_rows, ctx.rows_by_column[static_cast<size_t>(k)])) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    ctx.active[static_cast<size_t>(target)] = 0;
                    ++removed;
                }
                continue;
            }

            // budget >= 2: lower-bound check
            bool infeasible = false;
            double max_min_cost = 0.0;
            for (int r : uncovered_rows) {
                double cheapest = std::numeric_limits<double>::infinity();
                for (const CostCol &cc : columns_by_row[static_cast<size_t>(r)]) {
                    if (cc.col == target) continue;
                    if (!ctx.active[static_cast<size_t>(cc.col)]) continue;
                    if (cc.cost <= budget + tol) { cheapest = cc.cost; break; }
                }
                if (!std::isfinite(cheapest)) { infeasible = true; break; }
                if (cheapest > max_min_cost) max_min_cost = cheapest;
            }
            if (infeasible || max_min_cost > budget + tol) {
                ctx.active[static_cast<size_t>(target)] = 0;
                ++removed;
            }
        }
    done:
        return removed;
    }
};

} // anonymous namespace

std::vector<std::unique_ptr<IColumnPreprocessRule>>
make_preprocess_rules(const std::string &configured) {
    std::vector<std::unique_ptr<IColumnPreprocessRule>> rules;
    const auto tokens = splitCsvTokens(configured);
    if (tokens.empty()) {
        rules.push_back(std::make_unique<SingleColumnDominanceRule>());
        rules.push_back(std::make_unique<TwoColumnDominanceRule>());
        return rules;
    }
    for (const std::string &token : tokens) {
        if (token == "none") { rules.clear(); return rules; }
        if (token == "single_column_dominance" || token == "single")
            rules.push_back(std::make_unique<SingleColumnDominanceRule>());
        else if (token == "two_column_dominance" || token == "pair" || token == "two")
            rules.push_back(std::make_unique<TwoColumnDominanceRule>());
        else if (token == "cost_driven_replacement" || token == "cost_driven")
            rules.push_back(std::make_unique<CostDrivenReplacementRule>());
        else if (token == "incumbent_budget_pruning" || token == "incumbent_budget")
            rules.push_back(std::make_unique<IncumbentBudgetPruningRule>());
    }
    if (rules.empty()) {
        rules.push_back(std::make_unique<SingleColumnDominanceRule>());
        rules.push_back(std::make_unique<TwoColumnDominanceRule>());
    }
    return rules;
}

} // namespace scpsol
