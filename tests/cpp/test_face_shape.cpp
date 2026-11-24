// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <numbers>
#include <vector>

#include "origamizer/check.hpp"
#include "origamizer/mesh.hpp"
#include "support.hpp"

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;
using origamizer::face_shape;
using origamizer::PolygonSoup;
using origamizer::split_nonconvex_faces;
using test::flat_polygon;

namespace {

double piece_area(const origamizer::PolygonSoup& soup, std::size_t f) {
  return face_shape(soup, f).area;
}

}  // namespace

TEST_CASE("face_shape measures a unit square") {
  const auto soup = flat_polygon({{0, 0}, {1, 0}, {1, 1}, {0, 1}});
  const auto shape = face_shape(soup, 0);
  REQUIRE(shape.area == Approx(1.0));
  REQUIRE(shape.diameter == Approx(std::sqrt(2.0)));
  REQUIRE(shape.normal.z() == Approx(1.0));
  REQUIRE(shape.planarity_error == Approx(0.0).margin(1e-15));
  REQUIRE(shape.simple);
  REQUIRE(shape.convex);
}

TEST_CASE("face_shape orients the normal by the vertex order") {
  const auto soup = flat_polygon({{0, 0}, {0, 1}, {1, 1}, {1, 0}});
  REQUIRE(face_shape(soup, 0).normal.z() == Approx(-1.0));
  REQUIRE(face_shape(soup, 0).convex);
}

TEST_CASE("face_shape allows straight corners in a convex face") {
  const auto soup = flat_polygon({{0, 0}, {1, 0}, {2, 0}, {2, 2}, {0, 2}});
  REQUIRE(face_shape(soup, 0).convex);
}

TEST_CASE("face_shape measures nonplanarity relative to the diameter") {
  auto soup = flat_polygon({{0, 0}, {1, 0}, {1, 1}, {0, 1}});
  soup.points[2] = origamizer::Point3(1, 1, 0.2);
  const auto shape = face_shape(soup, 0);
  REQUIRE(shape.planarity_error > 1e-3);
  REQUIRE_THROWS_WITH(split_nonconvex_faces(soup), ContainsSubstring("face 0 is not planar"));
}

TEST_CASE("split_nonconvex_faces triangulates a nonplanar quadrilateral on request") {
  // A skew quad: lifting one corner makes it nonplanar. Diagonal 0-2 gives
  // triangles (0,1,2) and (0,2,3); diagonal 1-3 gives (0,1,3) and (1,2,3).
  auto soup = flat_polygon({{0, 0}, {4, 0}, {4, 1}, {0, 1}});
  soup.points[2] = origamizer::Point3(4, 1, 0.5);
  REQUIRE_THROWS_WITH(split_nonconvex_faces(soup), ContainsSubstring("set triangulate_nonplanar"));
  const auto split = split_nonconvex_faces(soup, {}, true);
  REQUIRE(split.triangulated == std::vector<std::size_t>{0});
  REQUIRE(split.split.empty());
  REQUIRE(split.soup.faces.size() == 2);
  REQUIRE(split.origin == std::vector<std::size_t>{0, 0});
  for (std::size_t f = 0; f < 2; ++f) {
    REQUIRE(split.soup.faces[f].size() == 3);
    REQUIRE(face_shape(split.soup, f).planarity_error < 1e-12);
    REQUIRE(face_shape(split.soup, f).normal.z() > 0);  // orientation kept
  }
}

TEST_CASE("triangulate_nonplanar picks the diagonal with the larger smallest angle") {
  // A kite-like skew quad where diagonal 1-3 avoids a sliver.
  PolygonSoup soup;
  soup.points = {origamizer::Point3(0, 0, 0), origamizer::Point3(1, -0.2, 0.05),
                 origamizer::Point3(2, 0, 0), origamizer::Point3(1, 3, 0.05)};
  soup.faces = {{0, 1, 2, 3}};
  REQUIRE(face_shape(soup, 0).planarity_error > 1e-6);
  const auto split = split_nonconvex_faces(soup, {}, true);
  REQUIRE(split.soup.faces[0] == std::vector<std::size_t>{0, 1, 3});
  REQUIRE(split.soup.faces[1] == std::vector<std::size_t>{1, 2, 3});
}

TEST_CASE("split_nonconvex_faces leaves convex faces alone") {
  const auto soup = flat_polygon({{0, 0}, {1, 0}, {1, 1}, {0, 1}});
  const auto split = split_nonconvex_faces(soup);
  REQUIRE(split.soup.faces == soup.faces);
  REQUIRE(split.origin == std::vector<std::size_t>{0});
  REQUIRE(split.split.empty());
}

TEST_CASE("split_nonconvex_faces cuts an L into two convex pieces along a diagonal") {
  const auto soup = flat_polygon({{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}});
  REQUIRE_FALSE(face_shape(soup, 0).convex);
  const auto split = split_nonconvex_faces(soup);
  REQUIRE(split.split == std::vector<std::size_t>{0});
  REQUIRE(split.soup.faces.size() == 2);
  REQUIRE(split.origin == std::vector<std::size_t>{0, 0});
  REQUIRE(split.soup.points == soup.points);  // diagonals only: no new vertices
  double total = 0;
  for (std::size_t f = 0; f < split.soup.faces.size(); ++f) {
    const auto shape = face_shape(split.soup, f);
    REQUIRE(shape.convex);
    REQUIRE(shape.normal.z() == Approx(1.0));  // pieces keep the orientation
    total += piece_area(split.soup, f);
  }
  REQUIRE(total == Approx(3.0));
}

TEST_CASE("split_nonconvex_faces handles a U and a straight corner next to a notch") {
  const auto u = flat_polygon({{0, 0}, {3, 0}, {3, 2}, {2, 2}, {2, 1}, {1, 1}, {1, 2}, {0, 2}});
  const auto split_u = split_nonconvex_faces(u);
  REQUIRE(split_u.soup.faces.size() == 3);
  const auto l = flat_polygon({{0, 0}, {1, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}});
  const auto split_l = split_nonconvex_faces(l);
  double total = 0;
  for (std::size_t f = 0; f < split_l.soup.faces.size(); ++f) {
    REQUIRE(face_shape(split_l.soup, f).convex);
    total += piece_area(split_l.soup, f);
  }
  REQUIRE(total == Approx(3.0));
}

TEST_CASE("split_nonconvex_faces rejects self-crossing and zero-area faces") {
  std::vector<std::pair<double, double>> star;
  for (int k = 0; k < 5; ++k) {
    const double t = 2 * std::numbers::pi * (2 * k) / 5;
    star.emplace_back(std::cos(t), std::sin(t));
  }
  REQUIRE_FALSE(face_shape(flat_polygon(star), 0).simple);
  REQUIRE_THROWS_WITH(split_nonconvex_faces(flat_polygon(star)),
                      ContainsSubstring("face 0 crosses itself"));
  REQUIRE_THROWS_WITH(split_nonconvex_faces(flat_polygon({{0, 0}, {1, 0}, {2, 0}})),
                      ContainsSubstring("zero area"));
}

TEST_CASE("prepare_mesh splits the nonconvex caps of the L prism") {
  const auto prepared =
      origamizer::prepare_mesh(origamizer::read_obj(test::example_path("l_prism")));
  REQUIRE(prepared.split_faces == std::vector<std::size_t>{0, 1});
  REQUIRE(prepared.mesh.number_of_faces() == 10);
  REQUIRE(prepared.face_origin.size() == 10);
  REQUIRE(origamizer::topology(prepared.mesh).euler_characteristic == 2);
  REQUIRE_FALSE(prepared.flipped_to_clean_side);
}
