// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <vector>

#include "origamizer/check.hpp"
#include "origamizer/mesh.hpp"
#include "support.hpp"

using Catch::Matchers::ContainsSubstring;
using origamizer::check_soup;
using origamizer::IndexPair;
using origamizer::InvariantError;
using origamizer::Point3;
using test::make_soup;

namespace {

origamizer::PolygonSoup tetrahedron() {
  return origamizer::read_obj(test::example_path("tetrahedron"));
}

/// Points in general position, so that only the combinatorics matter.
std::vector<Point3> scattered(std::size_t n) {
  std::vector<Point3> points;
  for (std::size_t i = 0; i < n; ++i) {
    const double t = static_cast<double>(i);
    points.emplace_back(std::cos(t), std::sin(t), 0.1 * t * t);
  }
  return points;
}

}  // namespace

TEST_CASE("check_soup accepts the example meshes") {
  for (const char* name : test::kExamples) {
    INFO(name);
    const auto report = check_soup(origamizer::read_obj(test::example_path(name)));
    REQUIRE(report.ok());
    REQUIRE(report.describe().empty());
    REQUIRE(report.components == 1);
  }
}

TEST_CASE("check_soup reports degenerate faces and unused vertices") {
  const auto report =
      check_soup(make_soup(scattered(5), {{0, 1, 2}, {0, 1, 1}, {0, 1}, {0, 2, 9}}));
  REQUIRE(report.degenerate_faces == std::vector<std::size_t>{1, 2, 3});
  REQUIRE(report.unused_vertices == std::vector<std::size_t>{3, 4});
  REQUIRE_FALSE(report.ok());
  REQUIRE_THAT(report.describe(), ContainsSubstring("3 degenerate faces"));
}

TEST_CASE("check_soup reports duplicate faces, including reversed copies") {
  const auto report = check_soup(make_soup(scattered(3), {{0, 1, 2}, {2, 0, 1}, {0, 2, 1}}));
  REQUIRE(report.duplicate_faces == std::vector<IndexPair>{{0, 1}, {0, 2}});
}

TEST_CASE("check_soup reports coincident vertices") {
  auto points = scattered(4);
  points[3] = points[0];
  const auto report = check_soup(make_soup(points, {{0, 1, 2}, {3, 2, 1}}));
  REQUIRE(report.coincident_vertices == std::vector<IndexPair>{{0, 3}});
}

TEST_CASE("check_soup reports an edge shared by three faces") {
  const auto report = check_soup(make_soup(scattered(5), {{0, 1, 2}, {1, 0, 3}, {0, 1, 4}}));
  REQUIRE(report.nonmanifold_edges == std::vector<IndexPair>{{0, 1}});
  REQUIRE(report.nonmanifold_vertices == std::vector<std::size_t>{0, 1});
}

TEST_CASE("check_soup reports a bowtie vertex and two components") {
  const auto report = check_soup(make_soup(scattered(5), {{0, 1, 2}, {0, 3, 4}}));
  REQUIRE(report.nonmanifold_vertices == std::vector<std::size_t>{0});
  REQUIRE(report.components == 2);
  REQUIRE_THAT(report.describe(), ContainsSubstring("2 connected components"));
}

TEST_CASE("check_soup detects an orientable soup with inconsistent faces and fixes it") {
  auto soup = make_soup(scattered(4), {{0, 1, 2}, {0, 1, 3}});
  const auto report = check_soup(soup);
  REQUIRE(report.ok_except_orientation());
  REQUIRE(report.orientable);
  REQUIRE_FALSE(report.consistently_oriented);
  REQUIRE(origamizer::make_orientation_consistent(soup) == std::vector<std::size_t>{1});
  REQUIRE(check_soup(soup).ok());
  REQUIRE(soup.faces[0] == std::vector<std::size_t>{0, 1, 2});
}

TEST_CASE("check_soup detects a nonorientable surface (5-vertex Moebius strip)") {
  auto soup = make_soup(scattered(5), {{0, 1, 2}, {1, 2, 3}, {2, 3, 4}, {3, 4, 0}, {4, 0, 1}});
  const auto report = check_soup(soup);
  REQUIRE(report.nonmanifold_edges.empty());
  REQUIRE(report.nonmanifold_vertices.empty());
  REQUIRE(report.components == 1);
  REQUIRE_FALSE(report.orientable);
  REQUIRE_THROWS_AS(origamizer::make_orientation_consistent(soup), InvariantError);
}

TEST_CASE("build_mesh keeps the vertex and face numbering") {
  const auto soup = tetrahedron();
  const auto mesh = origamizer::build_mesh(soup);
  REQUIRE(mesh.number_of_vertices() == 4);
  REQUIRE(mesh.number_of_faces() == 4);
  const auto back = origamizer::to_soup(mesh);
  REQUIRE(back.points == soup.points);
  for (std::size_t f = 0; f < 4; ++f) {
    auto a = back.faces[f], b = soup.faces[f];
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    REQUIRE(a == b);
  }
}

TEST_CASE("build_mesh rejects a soup that fails the checks and lists the problems") {
  REQUIRE_THROWS_WITH(origamizer::build_mesh(make_soup(scattered(5), {{0, 1, 2}, {0, 3, 4}})),
                      ContainsSubstring("not one connected") && ContainsSubstring("single fan"));
}
