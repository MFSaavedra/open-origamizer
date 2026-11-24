// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

#include "detail/quadratic_program.hpp"

#include <Highs.h>

#include <Eigen/OrderingMethods>
#include <Eigen/SparseCore>
#include <Eigen/SparseQR>
#include <algorithm>
#include <cmath>
#include <map>
#include <optional>
#include <utility>

#include "origamizer/check.hpp"

namespace origamizer::detail {

namespace {

/// A primal point with its row multipliers.
struct Point {
  std::vector<double> x;
  std::vector<double> y;
};

/// The objective (1/2) x^T D x + c^T x.
double objective_of(const QuadraticProgram& program, const std::vector<double>& x) {
  double value = 0;
  for (std::size_t j = 0; j < x.size(); ++j) {
    const double d = program.diagonal.empty() ? 0.0 : program.diagonal[j];
    value += (0.5 * d * x[j] + program.cost[j]) * x[j];
  }
  return value;
}

/// Checks a point against the KKT conditions of the program, independently of
/// HiGHS's status (docs/notes.tex, M2): primal feasibility, and with
/// z = D x + c - A^T y (stationarity defines z), the sign of every multiplier
/// y_i and z_j at its bound (dual feasibility and complementary slackness).
/// Returns what fails, or an empty string.
std::string kkt_failure(const QuadraticProgram& program, const HighsLp& lp, const Point& point) {
  constexpr double kRelative = 1e-6;
  const std::size_t n = program.variables();
  const std::size_t m = program.rows();
  const std::vector<double>& x = point.x;
  const std::vector<double>& y = point.y;

  std::vector<double> activity(m, 0.0), z(n, 0.0);
  double largest_gradient = 1.0;
  for (std::size_t j = 0; j < n; ++j) {
    z[j] = program.cost[j] + (program.diagonal.empty() ? 0.0 : program.diagonal[j] * x[j]);
    largest_gradient = std::max(largest_gradient, std::abs(z[j]));
    for (HighsInt k = lp.a_matrix_.start_[j]; k < lp.a_matrix_.start_[j + 1]; ++k) {
      const auto i = static_cast<std::size_t>(lp.a_matrix_.index_[static_cast<std::size_t>(k)]);
      const double a = lp.a_matrix_.value_[static_cast<std::size_t>(k)];
      activity[i] += a * x[j];
      z[j] -= a * y[i];
    }
  }
  const double dual_tolerance = kRelative * largest_gradient;

  // One bounded quantity (a variable or a row activity) with its multiplier.
  const auto check = [&](const char* what, std::size_t index, double value, double lo, double hi,
                         double multiplier) -> std::string {
    const double lo_tolerance = kRelative * (1.0 + std::abs(lo));
    const double hi_tolerance = kRelative * (1.0 + std::abs(hi));
    if (value < lo - lo_tolerance || value > hi + hi_tolerance) {
      return std::string(what) + " " + std::to_string(index) + " violates its bounds";
    }
    const bool at_lower = std::isfinite(lo) && value <= lo + lo_tolerance;
    const bool at_upper = std::isfinite(hi) && value >= hi - hi_tolerance;
    const bool sign_ok = (at_lower && at_upper) || (at_lower && multiplier >= -dual_tolerance) ||
                         (at_upper && multiplier <= dual_tolerance) ||
                         (!at_lower && !at_upper && std::abs(multiplier) <= dual_tolerance);
    if (!sign_ok) {
      return std::string(what) + " " + std::to_string(index) +
             " has a multiplier of the wrong sign (" + std::to_string(multiplier) + ")";
    }
    return {};
  };
  for (std::size_t j = 0; j < n; ++j) {
    if (auto failure = check("variable", j, x[j], program.lower[j], program.upper[j], z[j]);
        !failure.empty()) {
      return failure;
    }
  }
  for (std::size_t i = 0; i < m; ++i) {
    if (auto failure =
            check("row", i, activity[i], program.row_lower[i], program.row_upper[i], y[i]);
        !failure.empty()) {
      return failure;
    }
  }
  return {};
}

/// Polishing (docs/notes.tex, M2): keeps the active set of HiGHS's optimum
/// and solves its KKT equations to machine precision with a sparse QR:
///   D_FF x_F - A_RF^T y_R = -c_F,   A_RF x_F = b_R - A_RX x_X,
/// where X are the variables at an active bound (fixed there), F the others,
/// and R the active rows (equalities and rows at a bound) with targets b_R.
/// HiGHS's QP solver regularizes the Hessian and stops at tolerances near
/// 1e-7, which turns exact zeros and exact dependencies of the solution into
/// noise; polishing removes that noise. Returns nothing if the basis is
/// unusable or the equations cannot be solved; the caller still checks the
/// result.
std::optional<Point> polish(const QuadraticProgram& program, const HighsLp& lp,
                            const HighsBasis& basis) {
  const std::size_t n = program.variables();
  const std::size_t m = program.rows();
  if (!basis.valid || basis.col_status.size() != n || basis.row_status.size() != m) return {};

  constexpr std::size_t kNone = static_cast<std::size_t>(-1);
  Point point;
  point.x.assign(n, 0.0);
  point.y.assign(m, 0.0);
  std::vector<std::size_t> free_index(n, kNone);
  std::size_t free_count = 0;
  for (std::size_t j = 0; j < n; ++j) {
    switch (basis.col_status[j]) {
      case HighsBasisStatus::kLower:
        if (!std::isfinite(program.lower[j])) return {};
        point.x[j] = program.lower[j];
        break;
      case HighsBasisStatus::kUpper:
        if (!std::isfinite(program.upper[j])) return {};
        point.x[j] = program.upper[j];
        break;
      case HighsBasisStatus::kZero:  // a nonbasic free variable, at zero
        break;
      default:
        free_index[j] = free_count++;
    }
  }
  std::vector<std::size_t> active_index(m, kNone);
  std::vector<double> target;
  for (std::size_t i = 0; i < m; ++i) {
    double b = 0;
    if (program.row_lower[i] == program.row_upper[i] ||
        basis.row_status[i] == HighsBasisStatus::kLower) {
      b = program.row_lower[i];
    } else if (basis.row_status[i] == HighsBasisStatus::kUpper) {
      b = program.row_upper[i];
    } else {
      continue;
    }
    if (!std::isfinite(b)) return {};
    active_index[i] = target.size();
    target.push_back(b);
  }

  const auto size = static_cast<Eigen::Index>(free_count + target.size());
  if (size == 0) return point;
  std::vector<Eigen::Triplet<double>> triplets;
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(size);
  for (std::size_t j = 0; j < n; ++j) {
    const std::size_t f = free_index[j];
    if (f != kNone) {
      const double d = program.diagonal.empty() ? 0.0 : program.diagonal[j];
      if (d != 0) triplets.emplace_back(static_cast<int>(f), static_cast<int>(f), d);
      rhs(static_cast<Eigen::Index>(f)) = -program.cost[j];
    }
    for (HighsInt k = lp.a_matrix_.start_[j]; k < lp.a_matrix_.start_[j + 1]; ++k) {
      const auto i = static_cast<std::size_t>(lp.a_matrix_.index_[static_cast<std::size_t>(k)]);
      const double a = lp.a_matrix_.value_[static_cast<std::size_t>(k)];
      const std::size_t r = active_index[i];
      if (r == kNone) continue;
      const auto row = static_cast<int>(free_count + r);
      if (f != kNone) {
        triplets.emplace_back(row, static_cast<int>(f), a);
        triplets.emplace_back(static_cast<int>(f), row, -a);
      } else {
        target[r] -= a * point.x[j];
      }
    }
  }
  for (std::size_t r = 0; r < target.size(); ++r) {
    rhs(static_cast<Eigen::Index>(free_count + r)) = target[r];
  }
  Eigen::SparseMatrix<double> kkt(size, size);
  kkt.setFromTriplets(triplets.begin(), triplets.end());
  kkt.makeCompressed();
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> qr;
  qr.compute(kkt);
  if (qr.info() != Eigen::Success) return {};
  const Eigen::VectorXd solution = qr.solve(rhs);
  if (qr.info() != Eigen::Success || !solution.allFinite()) return {};
  // A singular system solved in the least-squares sense may miss: accept only
  // an exact solve.
  if ((kkt * solution - rhs).lpNorm<Eigen::Infinity>() >
      1e-9 * (1.0 + rhs.lpNorm<Eigen::Infinity>())) {
    return {};
  }
  for (std::size_t j = 0; j < n; ++j) {
    if (free_index[j] != kNone) point.x[j] = solution(static_cast<Eigen::Index>(free_index[j]));
  }
  for (std::size_t i = 0; i < m; ++i) {
    if (active_index[i] != kNone) {
      point.y[i] = solution(static_cast<Eigen::Index>(free_count + active_index[i]));
    }
  }
  return point;
}

}  // namespace

std::size_t QuadraticProgram::add_variable(double lo, double hi, double linear_cost) {
  cost.push_back(linear_cost);
  lower.push_back(lo);
  upper.push_back(hi);
  if (!diagonal.empty()) diagonal.push_back(0);
  return cost.size() - 1;
}

std::size_t QuadraticProgram::add_row(double lo, double hi) {
  row_lower.push_back(lo);
  row_upper.push_back(hi);
  return row_lower.size() - 1;
}

QpSolution solve(const QuadraticProgram& program) {
  const std::size_t n = program.variables();
  const std::size_t m = program.rows();
  ORIGAMIZER_CHECK(n > 0, "a program needs at least one variable");
  ORIGAMIZER_CHECK(program.lower.size() == n && program.upper.size() == n,
                   "variable bounds do not match the number of variables");
  ORIGAMIZER_CHECK(program.diagonal.empty() || program.diagonal.size() == n,
                   "the quadratic diagonal does not match the number of variables");
  ORIGAMIZER_CHECK(program.row_upper.size() == m, "row bounds do not match the number of rows");
  for (std::size_t k = 0; k < n; ++k) {
    ORIGAMIZER_CHECK(program.lower[k] <= program.upper[k], "a variable has lower > upper");
    ORIGAMIZER_CHECK(std::isfinite(program.cost[k]), "a cost is not finite");
    ORIGAMIZER_CHECK(program.diagonal.empty() || program.diagonal[k] >= 0,
                     "the quadratic part must be convex (nonnegative diagonal)");
  }

  // Column-wise sparse matrix; repeated entries are summed.
  std::map<std::pair<std::size_t, std::size_t>, double> summed;  // (column, row) -> value
  for (const auto& e : program.entries) {
    ORIGAMIZER_CHECK(e.row < m && e.column < n, "a matrix entry is out of range");
    ORIGAMIZER_CHECK(std::isfinite(e.value), "a matrix entry is not finite");
    summed[{e.column, e.row}] += e.value;
  }

  HighsModel model;
  HighsLp& lp = model.lp_;
  lp.num_col_ = static_cast<HighsInt>(n);
  lp.num_row_ = static_cast<HighsInt>(m);
  lp.col_cost_ = program.cost;
  lp.col_lower_ = program.lower;
  lp.col_upper_ = program.upper;
  lp.row_lower_ = program.row_lower;
  lp.row_upper_ = program.row_upper;
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = lp.num_col_;
  lp.a_matrix_.num_row_ = lp.num_row_;
  lp.a_matrix_.start_.assign(n + 1, 0);
  for (const auto& [key, value] : summed) {
    if (std::abs(value) <= kSmallMatrixValue) continue;
    lp.a_matrix_.index_.push_back(static_cast<HighsInt>(key.second));
    lp.a_matrix_.value_.push_back(value);
    ++lp.a_matrix_.start_[key.first + 1];
  }
  for (std::size_t k = 0; k < n; ++k) lp.a_matrix_.start_[k + 1] += lp.a_matrix_.start_[k];

  const bool quadratic =
      std::any_of(program.diagonal.begin(), program.diagonal.end(), [](double d) { return d > 0; });
  if (quadratic) {
    HighsHessian& hessian = model.hessian_;
    hessian.dim_ = lp.num_col_;
    hessian.format_ = HessianFormat::kTriangular;
    hessian.start_.assign(n + 1, 0);
    for (std::size_t k = 0; k < n; ++k) {
      if (program.diagonal[k] > 0) {
        hessian.index_.push_back(static_cast<HighsInt>(k));
        hessian.value_.push_back(program.diagonal[k]);
      }
      hessian.start_[k + 1] = static_cast<HighsInt>(hessian.index_.size());
    }
  }

  Highs highs;
  highs.setOptionValue("output_flag", false);
  highs.setOptionValue("threads", 1);
  QpSolution solution;
  if (const HighsStatus passed = highs.passModel(model); passed != HighsStatus::kOk) {
    solution.message = "HiGHS rejected the model (" + highsStatusToString(passed) + ")";
    return solution;
  }
  const HighsStatus run = highs.run();
  const HighsModelStatus status = highs.getModelStatus();
  solution.message = highs.modelStatusToString(status);
  if (run != HighsStatus::kError && status == HighsModelStatus::kOptimal) {
    const HighsSolution& found = highs.getSolution();
    if (!found.value_valid || !found.dual_valid) {
      solution.message = "HiGHS reported an optimum without duals";
      return solution;
    }
    // Prefer the polished point; fall back to HiGHS's own if polishing fails
    // or does not pass the KKT check. Either way, the point returned passes it.
    const Point unpolished{found.col_value, found.row_dual};
    const std::optional<Point> polished = polish(program, lp, highs.getBasis());
    std::string failure = "polishing failed";
    if (polished && (failure = kkt_failure(program, lp, *polished)).empty()) {
      solution.x = polished->x;
    } else if (const std::string own = kkt_failure(program, lp, unpolished); own.empty()) {
      solution.x = unpolished.x;
    } else {
      solution.message = "HiGHS reported an optimum that fails the KKT conditions: " + own +
                         " (polished: " + failure + ")";
      return solution;
    }
    solution.status = SolveStatus::kOptimal;
    solution.objective = objective_of(program, solution.x);
  } else if (status == HighsModelStatus::kInfeasible) {
    solution.status = SolveStatus::kInfeasible;
  } else if (status == HighsModelStatus::kUnbounded ||
             status == HighsModelStatus::kUnboundedOrInfeasible) {
    solution.status = SolveStatus::kUnbounded;
  }
  return solution;
}

}  // namespace origamizer::detail
