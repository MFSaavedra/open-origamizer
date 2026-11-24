// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

// Private helper: a sparse convex quadratic program and its solution with
// HiGHS. HiGHS types stay out of the public headers.

#pragma once

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace origamizer::detail {

inline constexpr double kInfinity = std::numeric_limits<double>::infinity();

/// HiGHS treats matrix entries of magnitude at most its small_matrix_value
/// (1e-9 by default) as zero and reports a warning when it drops them. The
/// wrapper drops them first, so that any warning from passModel still counts
/// as a rejection. In the mapping steps such entries are rounding errors,
/// such as cos(pi / 2) = 6e-17 in the translation rows.
inline constexpr double kSmallMatrixValue = 1e-9;

/// minimize   (1/2) sum_k diagonal[k] x_k^2 + cost^T x
/// subject to row_lower <= A x <= row_upper  and  lower <= x <= upper.
/// The quadratic part is diagonal with nonnegative entries (convex), which is
/// all the mapping steps need; an empty `diagonal` makes it a linear program.
struct QuadraticProgram {
  struct Entry {
    std::size_t row;
    std::size_t column;
    double value;
  };

  std::vector<double> cost;      ///< one per variable
  std::vector<double> diagonal;  ///< empty, or one nonnegative entry per variable
  std::vector<double> lower;     ///< one per variable; -kInfinity allowed
  std::vector<double> upper;     ///< one per variable; kInfinity allowed
  /// Nonzeros of A; repeated (row, column) pairs add up. After summing,
  /// entries of magnitude at most kSmallMatrixValue count as zero, so rows
  /// should be scaled to have coefficients of order one.
  std::vector<Entry> entries;
  std::vector<double> row_lower;
  std::vector<double> row_upper;

  std::size_t variables() const { return cost.size(); }
  std::size_t rows() const { return row_lower.size(); }
  /// Adds a variable with the given bounds and linear cost; returns its index.
  std::size_t add_variable(double lo, double hi, double linear_cost = 0);
  /// Adds an empty row with the given bounds; returns its index.
  std::size_t add_row(double lo, double hi);
};

enum class SolveStatus { kOptimal, kInfeasible, kUnbounded, kFailed };

struct QpSolution {
  SolveStatus status = SolveStatus::kFailed;
  std::vector<double> x;  ///< the solution when status is kOptimal
  double objective = 0;
  std::string message;  ///< the solver's model status, for error messages
};

/// Solves the program with HiGHS (output off, one thread, deterministic).
/// An optimum counts only if it passes an independent check of the KKT
/// conditions; otherwise the status is kFailed and the message says why.
/// Throws InvariantError if the input is malformed (sizes, indices, bounds).
QpSolution solve(const QuadraticProgram& program);

}  // namespace origamizer::detail
