// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

#include "origamizer/check.hpp"

#include <string>

namespace origamizer::detail {

void fail_check(std::string_view condition, std::string_view message,
                std::source_location where) {
  std::string text = "invariant failed: ";
  text += message;
  text += " [";
  text += condition;
  text += "] at ";
  text += where.file_name();
  text += ':';
  text += std::to_string(where.line());
  throw InvariantError(text);
}

}  // namespace origamizer::detail
