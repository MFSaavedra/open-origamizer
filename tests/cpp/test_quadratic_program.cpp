// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "detail/quadratic_program.hpp"
#include "origamizer/check.hpp"

using Catch::Approx;
using origamizer::detail::kInfinity;
using origamizer::detail::QuadraticProgram;
using origamizer::detail::solve;
using origamizer::detail::SolveStatus;

TEST_CASE("solve finds the optimum of a small linear program") {
  // minimize x + 2y subject to x + y = 1, x, y >= 0  ->  x = 1, y = 0
  QuadraticProgram lp;
  lp.add_variable(0, kInfinity, 1);
  lp.add_variable(0, kInfinity, 2);
  const auto row = lp.add_row(1, 1);
  lp.entries = {{row, 0, 1.0}, {row, 1, 1.0}};
  const auto s = solve(lp);
  REQUIRE(s.status == SolveStatus::kOptimal);
  REQUIRE(s.x[0] == Approx(1.0));
  REQUIRE(s.x[1] == Approx(0.0).margin(1e-12));
  REQUIRE(s.objective == Approx(1.0));
}

TEST_CASE("solve handles a diagonal quadratic objective") {
  // minimize (x^2 + y^2) / 2 subject to x + y = 1  ->  x = y = 1/2
  QuadraticProgram qp;
  qp.diagonal = {1.0, 1.0};
  qp.cost = {0.0, 0.0};
  qp.lower = {-kInfinity, -kInfinity};
  qp.upper = {kInfinity, kInfinity};
  const auto row = qp.add_row(1, 1);
  qp.entries = {{row, 0, 1.0}, {row, 1, 0.5}, {row, 1, 0.5}};  // repeated entries add up
  const auto s = solve(qp);
  REQUIRE(s.status == SolveStatus::kOptimal);
  REQUIRE(s.x[0] == Approx(0.5).epsilon(1e-6));
  REQUIRE(s.x[1] == Approx(0.5).epsilon(1e-6));
}

TEST_CASE("solve reports infeasible and unbounded programs") {
  QuadraticProgram infeasible;
  infeasible.add_variable(0, kInfinity, 1);
  const auto row = infeasible.add_row(-kInfinity, -1);  // x <= -1 contradicts x >= 0
  infeasible.entries = {{row, 0, 1.0}};
  REQUIRE(solve(infeasible).status == SolveStatus::kInfeasible);

  QuadraticProgram unbounded;
  unbounded.add_variable(-kInfinity, kInfinity, 1);  // minimize x with no lower bound
  REQUIRE(solve(unbounded).status == SolveStatus::kUnbounded);
}

TEST_CASE("solve rejects malformed programs") {
  QuadraticProgram bad;
  bad.add_variable(1, 0, 0);  // lower > upper
  REQUIRE_THROWS_AS(solve(bad), origamizer::InvariantError);
  QuadraticProgram out_of_range;
  out_of_range.add_variable(0, 1, 0);
  out_of_range.add_row(0, 1);
  out_of_range.entries = {{0, 5, 1.0}};
  REQUIRE_THROWS_AS(solve(out_of_range), origamizer::InvariantError);
}
