// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

/// \file
/// Always-on invariant checks.
///
/// Unlike assert(), ORIGAMIZER_CHECK stays active in release builds, wheels
/// included: a violated geometric invariant must stop the computation rather
/// than yield a silently wrong crease pattern (CLAUDE.md, Working style).

#pragma once

#include <source_location>
#include <stdexcept>
#include <string_view>

namespace origamizer {

/// Thrown when an invariant check fails. The Python bindings expose it as
/// origamizer.InvariantError, a subclass of RuntimeError.
class InvariantError : public std::logic_error {
 public:
  using std::logic_error::logic_error;
};

namespace detail {

/// Throws InvariantError with the message, the failed condition, and the
/// location of the check. Called only through ORIGAMIZER_CHECK.
[[noreturn]] void fail_check(std::string_view condition, std::string_view message,
                             std::source_location where);

}  // namespace detail
}  // namespace origamizer

/// Evaluates `condition` once, in every build type. If it is false, throws
/// origamizer::InvariantError carrying `message` (anything convertible to
/// std::string_view), the condition text, and the file and line of the check.
/// The macro is an expression, so it also works inside other expressions.
#define ORIGAMIZER_CHECK(condition, message)                         \
  (static_cast<bool>(condition)                                      \
       ? static_cast<void>(0)                                        \
       : ::origamizer::detail::fail_check(#condition, (message),     \
                                          std::source_location::current()))
