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

RowReductionResult row_reduce(
    int &nrows, int &ncols,
    std::vector<int> &csr_inds, std::vector<int> &csr_offs, std::vector<double> &csr_vals,
    std::vector<double> &obj, std::vector<int> &active_to_input,
    double tol, double time_limit_sec, int verbosity) {
    RowReductionResult result;
    if (nrows <= 0 || ncols <= 0) return result;

    auto deadline = std::chrono::steady_clock::time_point::max();
    if (time_limit_sec > 0.0) {
        deadline = std::chrono::steady_clock::now() +
                   std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                       std::chrono::duration<double>(time_limit_sec));
    }

    // Build adjacency
    std::vector<std::vector<int>> cols_by_row(static_cast<size_t>(nrows));
    std::vector<std::vector<int>> rows_by_col(static_cast<size_t>(ncols));
    for (int i = 0; i < nrows; ++i) {
        for (int k = csr_offs[static_cast<size_t>(i)];
             k < csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int col = csr_inds[static_cast<size_t>(k)];
            if (col >= 0 && col < ncols && csr_vals[static_cast<size_t>(k)] > tol) {
                cols_by_row[static_cast<size_t>(i)].push_back(col);
                rows_by_col[static_cast<size_t>(col)].push_back(i);
            }
        }
    }

    std::vector<char> col_active(static_cast<size_t>(ncols), 1);
    std::vector<char> row_active(static_cast<size_t>(nrows), 1);
    std::vector<char> col_fixed(static_cast<size_t>(ncols), 0);

    // ---- Phase 1: Essential column fixing (cascade) ----
    // A row covered by exactly 1 active column forces that column to 1.
    // Fixing a column removes all rows it covers, potentially creating new essentials.
    bool changed = true;
    while (changed) {
        changed = false;
        for (int i = 0; i < nrows; ++i) {
            if (!row_active[static_cast<size_t>(i)]) continue;
            int count = 0, last = -1;
            for (int c : cols_by_row[static_cast<size_t>(i)]) {
                if (col_active[static_cast<size_t>(c)]) {
                    ++count; last = c;
                    if (count > 1) break;
                }
            }
            if (count == 1) {
                col_fixed[static_cast<size_t>(last)] = 1;
                col_active[static_cast<size_t>(last)] = 0;
                result.fixed_cost += obj[static_cast<size_t>(last)];
                result.fixed_original_cols.push_back(
                    active_to_input[static_cast<size_t>(last)]);
                ++result.cols_fixed;
                for (int r : rows_by_col[static_cast<size_t>(last)]) {
                    if (row_active[static_cast<size_t>(r)]) {
                        row_active[static_cast<size_t>(r)] = 0;
                        ++result.rows_removed;
                        changed = true;
                    }
                }
            }
        }
    }

    // ---- Phase 2: Row domination ----
    // If the active columns covering row A are a subset of those covering row B,
    // then any solution satisfying A automatically satisfies B. Remove B.
    // Sort rows by coverage size (ascending) so smaller sets are checked first.
    std::vector<std::vector<int>> active_cols_sorted(static_cast<size_t>(nrows));
    std::vector<int> row_order;
    for (int i = 0; i < nrows; ++i) {
        if (!row_active[static_cast<size_t>(i)]) continue;
        for (int c : cols_by_row[static_cast<size_t>(i)]) {
            if (col_active[static_cast<size_t>(c)])
                active_cols_sorted[static_cast<size_t>(i)].push_back(c);
        }
        std::sort(active_cols_sorted[static_cast<size_t>(i)].begin(),
                  active_cols_sorted[static_cast<size_t>(i)].end());
        row_order.push_back(i);
    }
    std::sort(row_order.begin(), row_order.end(), [&](int a, int b) {
        return active_cols_sorted[static_cast<size_t>(a)].size() <
               active_cols_sorted[static_cast<size_t>(b)].size();
    });

    for (size_t ia = 0; ia < row_order.size(); ++ia) {
        if (std::chrono::steady_clock::now() >= deadline) break;
        const int a = row_order[ia];
        if (!row_active[static_cast<size_t>(a)]) continue;
        const auto &cols_a = active_cols_sorted[static_cast<size_t>(a)];
        for (size_t ib = ia + 1; ib < row_order.size(); ++ib) {
            const int b = row_order[ib];
            if (!row_active[static_cast<size_t>(b)]) continue;
            const auto &cols_b = active_cols_sorted[static_cast<size_t>(b)];
            if (cols_a.size() > cols_b.size()) continue;
            if (isSubsetSorted(cols_a, cols_b)) {
                row_active[static_cast<size_t>(b)] = 0;
                ++result.rows_removed;
            }
        }
    }

    // ---- Phase 3: Probing ----
    // For each active column j, simulate x_j = 0 and propagate forced fixings.
    // If the cascade leads to an uncoverable row, x_j must be 1.
    for (int j = 0; j < ncols; ++j) {
        if (std::chrono::steady_clock::now() >= deadline) break;
        if (!col_active[static_cast<size_t>(j)]) continue;

        // Skip columns that don't cover any active row
        bool covers_active = false;
        for (int r : rows_by_col[static_cast<size_t>(j)]) {
            if (row_active[static_cast<size_t>(r)]) { covers_active = true; break; }
        }
        if (!covers_active) continue;

        // Simulate x_j = 0 using rollback stack
        struct Change { int idx; bool is_col; };
        std::vector<Change> changes;
        col_active[static_cast<size_t>(j)] = 0;
        changes.push_back({j, true});

        // Worklist: columns just deactivated whose covered rows need checking
        std::vector<int> worklist;
        worklist.push_back(j);
        bool infeasible = false;

        while (!worklist.empty() && !infeasible) {
            const int deactivated = worklist.back();
            worklist.pop_back();
            for (int r : rows_by_col[static_cast<size_t>(deactivated)]) {
                if (!row_active[static_cast<size_t>(r)]) continue;
                int cnt = 0, last_c = -1;
                for (int c : cols_by_row[static_cast<size_t>(r)]) {
                    if (col_active[static_cast<size_t>(c)]) {
                        ++cnt; last_c = c;
                        if (cnt > 1) break;
                    }
                }
                if (cnt == 0) { infeasible = true; break; }
                if (cnt == 1) {
                    // Force last_c to 1: remove its rows, deactivate it
                    col_active[static_cast<size_t>(last_c)] = 0;
                    changes.push_back({last_c, true});
                    worklist.push_back(last_c);
                    for (int rr : rows_by_col[static_cast<size_t>(last_c)]) {
                        if (row_active[static_cast<size_t>(rr)]) {
                            row_active[static_cast<size_t>(rr)] = 0;
                            changes.push_back({rr, false});
                        }
                    }
                }
            }
        }

        // Rollback all temporary changes
        for (auto it = changes.rbegin(); it != changes.rend(); ++it) {
            if (it->is_col) col_active[static_cast<size_t>(it->idx)] = 1;
            else row_active[static_cast<size_t>(it->idx)] = 1;
        }

        if (infeasible) {
            // x_j = 0 is infeasible → fix x_j = 1
            col_fixed[static_cast<size_t>(j)] = 1;
            col_active[static_cast<size_t>(j)] = 0;
            result.fixed_cost += obj[static_cast<size_t>(j)];
            result.fixed_original_cols.push_back(
                active_to_input[static_cast<size_t>(j)]);
            ++result.cols_fixed;
            for (int r : rows_by_col[static_cast<size_t>(j)]) {
                if (row_active[static_cast<size_t>(r)]) {
                    row_active[static_cast<size_t>(r)] = 0;
                    ++result.rows_removed;
                }
            }
        }
    }

    // Re-run essential column cascade after probing fixings
    changed = true;
    while (changed) {
        changed = false;
        for (int i = 0; i < nrows; ++i) {
            if (!row_active[static_cast<size_t>(i)]) continue;
            int count = 0, last = -1;
            for (int c : cols_by_row[static_cast<size_t>(i)]) {
                if (col_active[static_cast<size_t>(c)]) {
                    ++count; last = c;
                    if (count > 1) break;
                }
            }
            if (count == 1) {
                col_fixed[static_cast<size_t>(last)] = 1;
                col_active[static_cast<size_t>(last)] = 0;
                result.fixed_cost += obj[static_cast<size_t>(last)];
                result.fixed_original_cols.push_back(
                    active_to_input[static_cast<size_t>(last)]);
                ++result.cols_fixed;
                for (int r : rows_by_col[static_cast<size_t>(last)]) {
                    if (row_active[static_cast<size_t>(r)]) {
                        row_active[static_cast<size_t>(r)] = 0;
                        ++result.rows_removed;
                        changed = true;
                    }
                }
            }
        }
    }

    // Remove columns that only cover inactive rows (dead columns)
    int dead_cols = 0;
    for (int j2 = 0; j2 < ncols; ++j2) {
        if (!col_active[static_cast<size_t>(j2)]) continue;
        bool any = false;
        for (int r : rows_by_col[static_cast<size_t>(j2)]) {
            if (row_active[static_cast<size_t>(r)]) { any = true; break; }
        }
        if (!any) { col_active[static_cast<size_t>(j2)] = 0; ++dead_cols; }
    }

    if (result.rows_removed == 0 && result.cols_fixed == 0 && dead_cols == 0)
        return result;

    // Rebuild CSR with active rows and active columns
    std::vector<int> old_col_to_new(static_cast<size_t>(ncols), -1);
    std::vector<int> new_active;
    int new_ncols = 0;
    for (int j2 = 0; j2 < ncols; ++j2) {
        if (col_active[static_cast<size_t>(j2)]) {
            old_col_to_new[static_cast<size_t>(j2)] = new_ncols;
            new_active.push_back(active_to_input[static_cast<size_t>(j2)]);
            ++new_ncols;
        }
    }

    std::vector<double> new_obj(static_cast<size_t>(new_ncols));
    for (int j2 = 0; j2 < ncols; ++j2) {
        const int nj = old_col_to_new[static_cast<size_t>(j2)];
        if (nj >= 0) new_obj[static_cast<size_t>(nj)] = obj[static_cast<size_t>(j2)];
    }

    std::vector<int> new_inds, new_offs;
    std::vector<double> new_vals;
    new_offs.push_back(0);
    int new_nrows = 0;
    for (int i = 0; i < nrows; ++i) {
        if (!row_active[static_cast<size_t>(i)]) continue;
        for (int k = csr_offs[static_cast<size_t>(i)];
             k < csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int c = csr_inds[static_cast<size_t>(k)];
            if (c >= 0 && c < ncols) {
                const int nc = old_col_to_new[static_cast<size_t>(c)];
                if (nc >= 0) {
                    new_inds.push_back(nc);
                    new_vals.push_back(csr_vals[static_cast<size_t>(k)]);
                }
            }
        }
        new_offs.push_back(static_cast<int>(new_vals.size()));
        ++new_nrows;
    }

    if (verbosity >= 3) {
        fprintf(stderr, "  Row reduction: %d rows, %d cols -> %d rows, %d cols "
                        "(%d essential, %d dead, cost %.6g)\n",
                nrows, ncols, new_nrows, new_ncols,
                result.cols_fixed, dead_cols, result.fixed_cost);
    }

    nrows = new_nrows;
    ncols = new_ncols;
    csr_inds = std::move(new_inds);
    csr_offs = std::move(new_offs);
    csr_vals = std::move(new_vals);
    obj = std::move(new_obj);
    active_to_input = std::move(new_active);

    return result;
}

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
