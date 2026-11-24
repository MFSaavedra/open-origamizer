// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

#include "detail/quadratic_program.hpp"

#include <Highs.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <utility>

#include "origamizer/check.hpp"

namespace origamizer::detail {

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
    if (value == 0) continue;
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
  if (highs.passModel(model) != HighsStatus::kOk) {
    solution.message = "HiGHS rejected the model";
    return solution;
  }
  const HighsStatus run = highs.run();
  const HighsModelStatus status = highs.getModelStatus();
  solution.message = highs.modelStatusToString(status);
  if (run != HighsStatus::kError && status == HighsModelStatus::kOptimal) {
    solution.status = SolveStatus::kOptimal;
    solution.x = highs.getSolution().col_value;
    solution.objective = highs.getInfo().objective_function_value;
  } else if (status == HighsModelStatus::kInfeasible) {
    solution.status = SolveStatus::kInfeasible;
  } else if (status == HighsModelStatus::kUnbounded ||
             status == HighsModelStatus::kUnboundedOrInfeasible) {
    solution.status = SolveStatus::kUnbounded;
  }
  return solution;
}

}  // namespace origamizer::detail
