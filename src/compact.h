#ifndef SCPSOL_COMPACT_H
#define SCPSOL_COMPACT_H

#include <algorithm>
#include <vector>

namespace scpsol {

// Build an old->new column map from keep flags (non-zero = keep).
// Returns the number of kept columns; dropped columns map to -1.
int build_column_map(const std::vector<char> &keep, std::vector<int> &old_to_new);

// Remove every column j with old_to_new[j] < 0 from a row-major CSR matrix and
// renumber the remaining entries. Entries whose column index is out of range
// are dropped as well.
void compact_csr_columns(int nrows, const std::vector<int> &old_to_new,
                         std::vector<int> &csr_inds, std::vector<int> &csr_offs,
                         std::vector<double> &csr_vals);

// Keep only the rows i with keep[i] != 0. Returns the new row count.
int compact_csr_rows(const std::vector<char> &keep,
                     std::vector<int> &csr_inds, std::vector<int> &csr_offs,
                     std::vector<double> &csr_vals);

// Apply an old->new column map to a per-column vector (costs, index maps, ...).
template <typename T>
void compact_column_vector(const std::vector<int> &old_to_new, int new_ncols,
                           std::vector<T> &values) {
    std::vector<T> out(static_cast<size_t>(new_ncols));
    const size_t n = std::min(values.size(), old_to_new.size());
    for (size_t j = 0; j < n; ++j) {
        const int nj = old_to_new[j];
        if (nj >= 0) out[static_cast<size_t>(nj)] = values[j];
    }
    values = std::move(out);
}

// The reduced set covering instance that root preprocessing operates on.
// Column j of the working model corresponds to column active_to_input[j]
// of the input instance.
struct WorkingModel {
    int nrows = 0;
    int ncols = 0;
    std::vector<int> csr_inds;
    std::vector<int> csr_offs;
    std::vector<double> csr_vals;
    std::vector<double> obj;
    std::vector<int> active_to_input;

    // Drop the columns with keep[j] == 0 and renumber the rest.
    // Returns the number of removed columns.
    int remove_columns(const std::vector<char> &keep);

    // Drop the rows with keep[i] == 0. Returns the number of removed rows.
    int remove_rows(const std::vector<char> &keep);
};

} // namespace scpsol

#endif // SCPSOL_COMPACT_H
