#include "decomposition.h"
#include "solver.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <numeric>
#include <queue>
#include <unordered_set>

namespace scpsol {

// ---- Connected components via BFS on bipartite row-column graph ----

// Returns component id per row (-1 for inactive rows).
// active_rows[i] = true means row i is in the residual problem.
// active_cols[j] = true means column j is in the residual problem.
static int find_components(
    const BaseRelaxationModel &base,
    const std::vector<bool> &active_rows,
    const std::vector<bool> &active_cols,
    std::vector<int> &row_comp) {

    const int nrows = base.nrows;
    const int ncols = base.ncols;
    row_comp.assign(static_cast<size_t>(nrows), -1);
    std::vector<int> col_comp(static_cast<size_t>(ncols), -1);

    // Build lightweight col->rows from CSR
    // (We avoid depending on cols_to_rows which may not be populated)
    std::vector<std::vector<int>> col_rows(static_cast<size_t>(ncols));
    for (int i = 0; i < nrows; ++i) {
        if (!active_rows[static_cast<size_t>(i)]) continue;
        for (int k = base.csr_offs[static_cast<size_t>(i)];
             k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int c = base.csr_inds[static_cast<size_t>(k)];
            if (c >= 0 && c < ncols && active_cols[static_cast<size_t>(c)])
                col_rows[static_cast<size_t>(c)].push_back(i);
        }
    }

    int num_components = 0;
    std::queue<int> q;

    for (int i = 0; i < nrows; ++i) {
        if (!active_rows[static_cast<size_t>(i)] ||
            row_comp[static_cast<size_t>(i)] >= 0)
            continue;

        // BFS from row i
        const int comp = num_components++;
        row_comp[static_cast<size_t>(i)] = comp;
        q.push(i);

        while (!q.empty()) {
            const int r = q.front(); q.pop();
            // Follow row r -> columns -> other rows
            for (int k = base.csr_offs[static_cast<size_t>(r)];
                 k < base.csr_offs[static_cast<size_t>(r) + 1]; ++k) {
                const int c = base.csr_inds[static_cast<size_t>(k)];
                if (c < 0 || c >= ncols || !active_cols[static_cast<size_t>(c)])
                    continue;
                if (col_comp[static_cast<size_t>(c)] >= 0) continue;
                col_comp[static_cast<size_t>(c)] = comp;
                for (int r2 : col_rows[static_cast<size_t>(c)]) {
                    if (row_comp[static_cast<size_t>(r2)] < 0) {
                        row_comp[static_cast<size_t>(r2)] = comp;
                        q.push(r2);
                    }
                }
            }
        }
    }
    return num_components;
}

// ---- Linking column heuristic ----

DecompositionState compute_decomposition(
    const BaseRelaxationModel &base,
    int max_linkers,
    int verbosity) {

    DecompositionState state;
    state.computed = true;

    const int nrows = base.nrows;
    const int ncols = base.ncols;
    if (nrows <= 0 || ncols <= 0) return state;

    // Start with all rows/cols active
    std::vector<bool> active_rows(static_cast<size_t>(nrows), true);
    std::vector<bool> active_cols(static_cast<size_t>(ncols), true);

    // Build row sets per column for scoring
    std::vector<std::vector<int>> col_rows(static_cast<size_t>(ncols));
    for (int i = 0; i < nrows; ++i) {
        for (int k = base.csr_offs[static_cast<size_t>(i)];
             k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int c = base.csr_inds[static_cast<size_t>(k)];
            if (c >= 0 && c < ncols)
                col_rows[static_cast<size_t>(c)].push_back(i);
        }
    }

    // Check initial connectivity
    std::vector<int> row_comp;
    int num_comp = find_components(base, active_rows, active_cols, row_comp);
    if (num_comp >= 2) {
        // Already decomposable without any linking columns
        state.enabled = true;
        // Build blocks
        state.blocks.resize(static_cast<size_t>(num_comp));
        for (int i = 0; i < nrows; ++i) {
            if (row_comp[static_cast<size_t>(i)] >= 0)
                state.blocks[static_cast<size_t>(row_comp[static_cast<size_t>(i)])].rows.push_back(i);
        }
        // Assign columns to blocks
        for (int j = 0; j < ncols; ++j) {
            int block = -1;
            for (int r : col_rows[static_cast<size_t>(j)]) {
                const int rc = row_comp[static_cast<size_t>(r)];
                if (rc >= 0) { block = rc; break; }
            }
            if (block >= 0)
                state.blocks[static_cast<size_t>(block)].cols.push_back(j);
        }
        if (verbosity >= 2)
            fprintf(stderr, "  Decomposition: %d blocks (no linking needed)\n", num_comp);
        return state;
    }

    // Greedy: find columns whose row-removal creates the most/smallest components
    std::vector<int> linking;
    std::unordered_set<int> linking_set;

    // Count initial active rows
    int initial_active_rows = 0;
    for (int i = 0; i < nrows; ++i)
        if (active_rows[static_cast<size_t>(i)]) ++initial_active_rows;

    for (int step = 0; step < max_linkers; ++step) {
        int best_col = -1;
        int best_num_comp = 0;
        int best_largest = nrows + 1;
        int best_active_rows_left = nrows + 1; // for tiebreaking when nc=1

        for (int j = 0; j < ncols; ++j) {
            if (linking_set.count(j)) continue;
            // Count active rows this column covers
            int covers = 0;
            for (int r : col_rows[static_cast<size_t>(j)])
                if (active_rows[static_cast<size_t>(r)]) ++covers;
            if (covers == 0) continue;

            // Simulate: deactivate column j and its rows
            std::vector<bool> test_rows = active_rows;
            std::vector<bool> test_cols = active_cols;
            test_cols[static_cast<size_t>(j)] = false;
            for (int r : col_rows[static_cast<size_t>(j)]) {
                if (active_rows[static_cast<size_t>(r)])
                    test_rows[static_cast<size_t>(r)] = false;
            }

            std::vector<int> test_comp;
            int nc = find_components(base, test_rows, test_cols, test_comp);

            // Compute largest component and active rows left
            int largest = 0;
            int rows_left = 0;
            if (nc >= 2) {
                std::vector<int> comp_size(static_cast<size_t>(nc), 0);
                for (int i = 0; i < nrows; ++i)
                    if (test_comp[static_cast<size_t>(i)] >= 0)
                        ++comp_size[static_cast<size_t>(test_comp[static_cast<size_t>(i)])];
                for (int s : comp_size) { largest = std::max(largest, s); rows_left += s; }
            } else {
                for (int i = 0; i < nrows; ++i)
                    if (test_comp[static_cast<size_t>(i)] >= 0) ++rows_left;
                largest = rows_left;
            }

            // Score: prefer nc >= 2 (splits); among nc=1, prefer most rows removed
            if (nc > best_num_comp ||
                (nc == best_num_comp && nc >= 2 && largest < best_largest) ||
                (nc == best_num_comp && nc == 1 && rows_left < best_active_rows_left)) {
                best_col = j;
                best_num_comp = nc;
                best_largest = largest;
                best_active_rows_left = rows_left;
            }
        }

        if (best_col < 0) break; // no column creates a split

        linking.push_back(best_col);
        linking_set.insert(best_col);
        // Deactivate this column and its rows
        active_cols[static_cast<size_t>(best_col)] = false;
        for (int r : col_rows[static_cast<size_t>(best_col)])
            active_rows[static_cast<size_t>(r)] = false;

        // Check if we have >=2 components now
        num_comp = find_components(base, active_rows, active_cols, row_comp);
        if (verbosity >= 3)
            fprintf(stderr, "  Decomp step %d: link col %d -> %d components\n",
                    step + 1, best_col, num_comp);
        if (num_comp >= 2) break;
    }

    if (num_comp < 2) {
        if (verbosity >= 2)
            fprintf(stderr, "  Decomposition: no split found (%d linkers tried)\n",
                    static_cast<int>(linking.size()));
        return state; // enabled stays false
    }

    state.enabled = true;
    state.linking_cols = linking;

    // Build blocks from final component assignment
    state.blocks.resize(static_cast<size_t>(num_comp));
    for (int i = 0; i < nrows; ++i) {
        if (row_comp[static_cast<size_t>(i)] >= 0)
            state.blocks[static_cast<size_t>(row_comp[static_cast<size_t>(i)])].rows.push_back(i);
    }
    for (int j = 0; j < ncols; ++j) {
        if (!active_cols[static_cast<size_t>(j)]) continue; // linking col
        int block = -1;
        for (int r : col_rows[static_cast<size_t>(j)]) {
            const int rc = row_comp[static_cast<size_t>(r)];
            if (rc >= 0) { block = rc; break; }
        }
        if (block >= 0)
            state.blocks[static_cast<size_t>(block)].cols.push_back(j);
    }

    if (verbosity >= 2) {
        fprintf(stderr, "  Decomposition: %d blocks, %d linking cols [",
                num_comp, static_cast<int>(linking.size()));
        for (size_t i = 0; i < linking.size(); ++i)
            fprintf(stderr, "%s%d", i ? "," : "", linking[i]);
        fprintf(stderr, "]\n");
        for (int b = 0; b < num_comp; ++b)
            fprintf(stderr, "    Block %d: %d rows, %d cols\n",
                    b, static_cast<int>(state.blocks[static_cast<size_t>(b)].rows.size()),
                    static_cast<int>(state.blocks[static_cast<size_t>(b)].cols.size()));
    }

    return state;
}

// ---- Query helpers ----

bool all_linkers_fixed(
    const DecompositionState &state,
    const std::vector<BranchDecision> &decisions) {

    if (state.linking_cols.empty()) return true;

    std::unordered_set<int> fixed_vars;
    for (const auto &d : decisions)
        fixed_vars.insert(d.var_index);

    for (int lc : state.linking_cols) {
        if (fixed_vars.find(lc) == fixed_vars.end())
            return false;
    }
    return true;
}

int pick_linking_branch_var(
    const DecompositionState &state,
    const std::vector<BranchDecision> &decisions,
    const std::vector<int> &fractional) {

    std::unordered_set<int> fixed_vars;
    for (const auto &d : decisions)
        fixed_vars.insert(d.var_index);

    std::unordered_set<int> frac_set(fractional.begin(), fractional.end());

    for (int lc : state.linking_cols) {
        if (fixed_vars.find(lc) == fixed_vars.end() &&
            frac_set.find(lc) != frac_set.end())
            return lc;
    }
    return -1;
}

// ---- Decomposed solve ----

DecompositionResult solve_decomposed(
    const BaseRelaxationModel &base,
    const std::vector<BranchDecision> &decisions,
    const DecompositionState &state,
    double incumbent_obj,
    const SolverConfig &config,
    double time_remaining,
    int verbosity) {

    DecompositionResult result;
    const int ncols = base.ncols;
    const int nrows = base.nrows;

    // Determine which rows are covered by linking cols fixed to 1
    std::vector<int8_t> col_fix(static_cast<size_t>(ncols), 0);
    for (const auto &d : decisions) {
        if (d.var_index >= 0 && d.var_index < ncols)
            col_fix[static_cast<size_t>(d.var_index)] =
                static_cast<int8_t>(d.fix_value == 1 ? 1 : -1);
    }

    // Cost of linking columns fixed to 1
    double linking_cost = 0.0;
    for (int lc : state.linking_cols) {
        if (col_fix[static_cast<size_t>(lc)] == 1)
            linking_cost += base.obj[static_cast<size_t>(lc)];
    }

    // Rebuild active rows: rows not covered by any col fixed to 1
    // (for correctness, re-derive blocks at this node instead of using
    //  precomputed blocks, since propagation may have fixed more columns)
    std::vector<bool> row_covered(static_cast<size_t>(nrows), false);
    for (int j = 0; j < ncols; ++j) {
        if (col_fix[static_cast<size_t>(j)] != 1) continue;
        for (int k = base.csr_offs[0]; k < static_cast<int>(base.csr_inds.size()); ++k) {
            // Need to find rows covered by column j via CSR
        }
    }

    // Use cols_to_rows if available, otherwise scan CSR
    // Mark rows covered by any column fixed to 1
    for (int i = 0; i < nrows; ++i) {
        for (int k = base.csr_offs[static_cast<size_t>(i)];
             k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int c = base.csr_inds[static_cast<size_t>(k)];
            if (c >= 0 && c < ncols && col_fix[static_cast<size_t>(c)] == 1)
                row_covered[static_cast<size_t>(i)] = true;
        }
    }

    // Build active rows/cols for component detection
    std::vector<bool> active_rows(static_cast<size_t>(nrows));
    std::vector<bool> active_cols(static_cast<size_t>(ncols));
    for (int i = 0; i < nrows; ++i)
        active_rows[static_cast<size_t>(i)] = !row_covered[static_cast<size_t>(i)];
    for (int j = 0; j < ncols; ++j)
        active_cols[static_cast<size_t>(j)] = (col_fix[static_cast<size_t>(j)] == 0);

    // Find connected components of residual problem
    std::vector<int> row_comp;
    int num_comp = find_components(base, active_rows, active_cols, row_comp);

    if (num_comp < 2) {
        // No decomposition possible at this node
        return result;
    }

    // Build column-to-block mapping
    std::vector<int> col_block(static_cast<size_t>(ncols), -1);
    // Build col_rows for active columns
    std::vector<std::vector<int>> col_rows(static_cast<size_t>(ncols));
    for (int i = 0; i < nrows; ++i) {
        if (!active_rows[static_cast<size_t>(i)]) continue;
        for (int k = base.csr_offs[static_cast<size_t>(i)];
             k < base.csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int c = base.csr_inds[static_cast<size_t>(k)];
            if (c >= 0 && c < ncols && active_cols[static_cast<size_t>(c)])
                col_rows[static_cast<size_t>(c)].push_back(i);
        }
    }
    for (int j = 0; j < ncols; ++j) {
        if (!active_cols[static_cast<size_t>(j)]) continue;
        for (int r : col_rows[static_cast<size_t>(j)]) {
            if (row_comp[static_cast<size_t>(r)] >= 0) {
                col_block[static_cast<size_t>(j)] = row_comp[static_cast<size_t>(r)];
                break;
            }
        }
    }

    if (verbosity >= 2)
        fprintf(stderr, "  Decomposed solve: %d blocks, linking cost=%.0f\n",
                num_comp, linking_cost);

    // Solve each block as an independent sub-problem
    double total_sub_obj = linking_cost;
    result.combined_solution.assign(static_cast<size_t>(ncols), 0.0);

    // Mark linking cols fixed to 1 in combined solution
    for (int lc : state.linking_cols) {
        if (col_fix[static_cast<size_t>(lc)] == 1)
            result.combined_solution[static_cast<size_t>(lc)] = 1.0;
    }
    // Also mark any non-linking col fixed to 1
    for (const auto &d : decisions) {
        if (d.fix_value == 1 && d.var_index >= 0 && d.var_index < ncols)
            result.combined_solution[static_cast<size_t>(d.var_index)] = 1.0;
    }

    double time_per_block = time_remaining / std::max(1, num_comp);

    for (int b = 0; b < num_comp; ++b) {
        // Collect rows and columns for this block
        std::vector<int> block_rows, block_cols;
        for (int i = 0; i < nrows; ++i) {
            if (row_comp[static_cast<size_t>(i)] == b)
                block_rows.push_back(i);
        }
        for (int j = 0; j < ncols; ++j) {
            if (col_block[static_cast<size_t>(j)] == b)
                block_cols.push_back(j);
        }

        if (block_rows.empty()) continue;

        if (block_cols.empty()) {
            // No columns to cover these rows — infeasible block
            if (verbosity >= 2)
                fprintf(stderr, "    Block %d: infeasible (no columns)\n", b);
            return result; // solved = false
        }

        // Build ScpInstance for this block
        ScpInstance sub;
        sub.nrows = static_cast<int>(block_rows.size());
        sub.ncols = static_cast<int>(block_cols.size());

        // Column index mapping: block_col_idx -> parent col
        std::vector<int> sub_to_parent(block_cols.begin(), block_cols.end());
        std::vector<int> parent_to_sub(static_cast<size_t>(ncols), -1);
        for (int jj = 0; jj < sub.ncols; ++jj)
            parent_to_sub[static_cast<size_t>(sub_to_parent[static_cast<size_t>(jj)])] = jj;

        // Costs
        sub.costs.resize(static_cast<size_t>(sub.ncols));
        for (int jj = 0; jj < sub.ncols; ++jj)
            sub.costs[static_cast<size_t>(jj)] =
                base.obj[static_cast<size_t>(sub_to_parent[static_cast<size_t>(jj)])];

        // Build CSR for sub-problem
        sub.csr_offsets.resize(static_cast<size_t>(sub.nrows) + 1, 0);
        for (int ii = 0; ii < sub.nrows; ++ii) {
            const int parent_row = block_rows[static_cast<size_t>(ii)];
            for (int k = base.csr_offs[static_cast<size_t>(parent_row)];
                 k < base.csr_offs[static_cast<size_t>(parent_row) + 1]; ++k) {
                const int pc = base.csr_inds[static_cast<size_t>(k)];
                if (pc >= 0 && pc < ncols &&
                    parent_to_sub[static_cast<size_t>(pc)] >= 0) {
                    sub.csr_indices.push_back(parent_to_sub[static_cast<size_t>(pc)]);
                    sub.csr_values.push_back(base.csr_vals[static_cast<size_t>(k)]);
                }
            }
            sub.csr_offsets[static_cast<size_t>(ii) + 1] =
                static_cast<int>(sub.csr_indices.size());
        }

        // Configure sub-solver with reduced limits
        SolverConfig sub_config = config;
        sub_config.verbosity = std::max(0, verbosity - 1);
        sub_config.time_limit_seconds = time_per_block;
        sub_config.max_nodes = std::min(config.max_nodes, 50000);
        sub_config.decomposition_mode = "off"; // no recursive decomposition

        if (verbosity >= 2)
            fprintf(stderr, "    Block %d: %d rows, %d cols -> solving...\n",
                    b, sub.nrows, sub.ncols);

        SolverResult sub_result = solve(sub, sub_config);

        if (sub_result.status == "Infeasible" ||
            !std::isfinite(sub_result.primal_obj)) {
            if (verbosity >= 2)
                fprintf(stderr, "    Block %d: infeasible\n", b);
            return result; // solved = false
        }

        total_sub_obj += sub_result.primal_obj;

        // Map sub-problem solution back to parent columns
        for (int jj = 0; jj < sub.ncols; ++jj) {
            if (static_cast<size_t>(jj) < sub_result.solution.size() &&
                sub_result.solution[static_cast<size_t>(jj)] > 0.5) {
                result.combined_solution[static_cast<size_t>(
                    sub_to_parent[static_cast<size_t>(jj)])] = 1.0;
            }
        }

        if (verbosity >= 2)
            fprintf(stderr, "    Block %d: %s obj=%.0f (%d nodes)\n",
                    b, sub_result.status.c_str(), sub_result.primal_obj,
                    sub_result.nodes_processed);
    }

    result.solved = true;
    result.combined_obj = total_sub_obj;
    result.improved = (total_sub_obj < incumbent_obj - 1e-6);

    if (verbosity >= 2)
        fprintf(stderr, "  Decomposed total: %.0f (incumbent: %.0f)%s\n",
                total_sub_obj, incumbent_obj,
                result.improved ? " -> improved!" : "");

    return result;
}

} // namespace scpsol
