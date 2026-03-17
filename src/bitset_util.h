#ifndef SCPSOL_BITSET_UTIL_H
#define SCPSOL_BITSET_UTIL_H

#include <cstdint>
#include <vector>

namespace scpsol {

// Fixed-width bitset for fast subset/union operations in preprocessing.
// Stores bits in a vector of uint64_t words.
struct DynBitset {
    std::vector<uint64_t> words;
    int nwords = 0;

    void init(int nbits) {
        nwords = (nbits + 63) / 64;
        words.assign(static_cast<size_t>(nwords), 0ULL);
    }

    void set(int bit) {
        words[static_cast<size_t>(bit / 64)] |= (1ULL << (bit % 64));
    }

    // Returns true if this is a subset of other (every bit set in this is also set in other)
    bool is_subset_of(const DynBitset &other) const {
        for (int i = 0; i < nwords; ++i) {
            if ((words[static_cast<size_t>(i)] & other.words[static_cast<size_t>(i)]) !=
                words[static_cast<size_t>(i)])
                return false;
        }
        return true;
    }

    // Returns true if this is a subset of (a | b)
    bool is_subset_of_union(const DynBitset &a, const DynBitset &b) const {
        for (int i = 0; i < nwords; ++i) {
            uint64_t u = a.words[static_cast<size_t>(i)] | b.words[static_cast<size_t>(i)];
            if ((words[static_cast<size_t>(i)] & u) != words[static_cast<size_t>(i)])
                return false;
        }
        return true;
    }

    // Returns true if this is a subset of (a | b | c)
    bool is_subset_of_union3(const DynBitset &a, const DynBitset &b, const DynBitset &c) const {
        for (int i = 0; i < nwords; ++i) {
            uint64_t u = a.words[static_cast<size_t>(i)] |
                         b.words[static_cast<size_t>(i)] |
                         c.words[static_cast<size_t>(i)];
            if ((words[static_cast<size_t>(i)] & u) != words[static_cast<size_t>(i)])
                return false;
        }
        return true;
    }

    int popcount() const {
        int count = 0;
        for (int i = 0; i < nwords; ++i) {
#if defined(__GNUC__) || defined(__clang__)
            count += __builtin_popcountll(words[static_cast<size_t>(i)]);
#else
            uint64_t v = words[static_cast<size_t>(i)];
            while (v) { ++count; v &= v - 1; }
#endif
        }
        return count;
    }
};

// Build bitsets for all columns from rows_by_column.
inline std::vector<DynBitset> build_column_bitsets(
    int ncols, int nrows,
    const std::vector<std::vector<int>> &rows_by_column,
    const std::vector<char> &active) {
    std::vector<DynBitset> bitsets(static_cast<size_t>(ncols));
    for (int j = 0; j < ncols; ++j) {
        bitsets[static_cast<size_t>(j)].init(nrows);
        if (!active[static_cast<size_t>(j)]) continue;
        for (int row : rows_by_column[static_cast<size_t>(j)]) {
            bitsets[static_cast<size_t>(j)].set(row);
        }
    }
    return bitsets;
}

} // namespace scpsol

#endif // SCPSOL_BITSET_UTIL_H
