#ifndef SCPSOL_READER_H
#define SCPSOL_READER_H

#include "model.h"
#include <string>

namespace scpsol {

ScpInstance read_scp_file(const std::string &path);

} // namespace scpsol

#endif // SCPSOL_READER_H
