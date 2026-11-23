// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

// Confirms that the CGAL and Eigen dependencies are wired into the build.

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

TEST_CASE("CGAL kernel predicates are available") {
  const K::Point_2 a(0, 0), b(1, 0), c(0, 1);
  REQUIRE(CGAL::orientation(a, b, c) == CGAL::LEFT_TURN);
  REQUIRE(CGAL::orientation(a, b, K::Point_2(2, 0)) == CGAL::COLLINEAR);
}

TEST_CASE("Eigen dense algebra is available") {
  Eigen::Matrix2d m;
  m << 2, 1, 1, 3;
  REQUIRE(m.determinant() == 5.0);
}
