#include "reader.h"

#include <fstream>
#include <stdexcept>

namespace scpsol {

ScpInstance read_scp_file(const std::string &path) {
    ScpInstance inst;
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Failed to open input file: " + path);
    }

    in >> inst.nrows >> inst.ncols;
    if (!in || inst.nrows <= 0 || inst.ncols <= 0) {
        throw std::runtime_error("Failed to parse SCP model dimensions");
    }

    inst.costs.resize(static_cast<size_t>(inst.ncols));
    for (int j = 0; j < inst.ncols; ++j) {
        if (!(in >> inst.costs[static_cast<size_t>(j)])) {
            throw std::runtime_error("Failed to parse objective coefficient at column " + std::to_string(j));
        }
    }

    inst.csr_offsets.push_back(0);
    for (int i = 0; i < inst.nrows; ++i) {
        int count;
        if (!(in >> count)) {
            throw std::runtime_error("Failed to parse row " + std::to_string(i) + " column count");
        }
        for (int j = 0; j < count; ++j) {
            int idx;
            if (!(in >> idx)) {
                throw std::runtime_error("Failed to parse column index at row " + std::to_string(i));
            }
            inst.csr_indices.push_back(idx - 1); // 1-indexed to 0-indexed
            inst.csr_values.push_back(1.0);
        }
        inst.csr_offsets.push_back(static_cast<int>(inst.csr_values.size()));
    }

    return inst;
}

} // namespace scpsol
