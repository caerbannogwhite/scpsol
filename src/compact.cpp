#include "compact.h"

namespace scpsol {

int build_column_map(const std::vector<char> &keep, std::vector<int> &old_to_new) {
    old_to_new.assign(keep.size(), -1);
    int next = 0;
    for (size_t j = 0; j < keep.size(); ++j) {
        if (keep[j]) old_to_new[j] = next++;
    }
    return next;
}

void compact_csr_columns(int nrows, const std::vector<int> &old_to_new,
                         std::vector<int> &csr_inds, std::vector<int> &csr_offs,
                         std::vector<double> &csr_vals) {
    const int ncols_old = static_cast<int>(old_to_new.size());
    std::vector<int> new_inds;
    std::vector<double> new_vals;
    std::vector<int> new_offs;
    new_inds.reserve(csr_inds.size());
    new_vals.reserve(csr_vals.size());
    new_offs.reserve(static_cast<size_t>(nrows) + 1);
    new_offs.push_back(0);
    for (int i = 0; i < nrows; ++i) {
        for (int k = csr_offs[static_cast<size_t>(i)]; k < csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            const int c = csr_inds[static_cast<size_t>(k)];
            if (c < 0 || c >= ncols_old) continue;
            const int m = old_to_new[static_cast<size_t>(c)];
            if (m < 0) continue;
            new_inds.push_back(m);
            new_vals.push_back(csr_vals[static_cast<size_t>(k)]);
        }
        new_offs.push_back(static_cast<int>(new_vals.size()));
    }
    csr_inds = std::move(new_inds);
    csr_offs = std::move(new_offs);
    csr_vals = std::move(new_vals);
}

int compact_csr_rows(const std::vector<char> &keep,
                     std::vector<int> &csr_inds, std::vector<int> &csr_offs,
                     std::vector<double> &csr_vals) {
    const int nrows = static_cast<int>(keep.size());
    std::vector<int> new_inds;
    std::vector<double> new_vals;
    std::vector<int> new_offs;
    new_inds.reserve(csr_inds.size());
    new_vals.reserve(csr_vals.size());
    new_offs.reserve(static_cast<size_t>(nrows) + 1);
    new_offs.push_back(0);
    int new_nrows = 0;
    for (int i = 0; i < nrows; ++i) {
        if (!keep[static_cast<size_t>(i)]) continue;
        for (int k = csr_offs[static_cast<size_t>(i)]; k < csr_offs[static_cast<size_t>(i) + 1]; ++k) {
            new_inds.push_back(csr_inds[static_cast<size_t>(k)]);
            new_vals.push_back(csr_vals[static_cast<size_t>(k)]);
        }
        new_offs.push_back(static_cast<int>(new_vals.size()));
        ++new_nrows;
    }
    csr_inds = std::move(new_inds);
    csr_offs = std::move(new_offs);
    csr_vals = std::move(new_vals);
    return new_nrows;
}

int WorkingModel::remove_columns(const std::vector<char> &keep) {
    std::vector<int> old_to_new;
    const int new_ncols = build_column_map(keep, old_to_new);
    const int removed = ncols - new_ncols;
    if (removed <= 0) return 0;
    compact_csr_columns(nrows, old_to_new, csr_inds, csr_offs, csr_vals);
    compact_column_vector(old_to_new, new_ncols, obj);
    compact_column_vector(old_to_new, new_ncols, active_to_input);
    ncols = new_ncols;
    return removed;
}

int WorkingModel::remove_rows(const std::vector<char> &keep) {
    const int new_nrows = compact_csr_rows(keep, csr_inds, csr_offs, csr_vals);
    const int removed = nrows - new_nrows;
    nrows = new_nrows;
    return removed;
}

} // namespace scpsol
