// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <numbers>
#include <string>

#include "origamizer/check.hpp"
#include "origamizer/mesh.hpp"
#include "support.hpp"

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;
using origamizer::gaussian_curvature;
using origamizer::SurfaceMesh;
using origamizer::topology;
using std::numbers::pi;

namespace {

SurfaceMesh example(const std::string& name) {
  return origamizer::prepare_mesh(origamizer::read_obj(test::example_path(name))).mesh;
}

double total_curvature(const SurfaceMesh& mesh) {
  double sum = 0;
  for (const auto v : mesh.vertices()) sum += gaussian_curvature(mesh, v);
  return sum;
}

}  // namespace

TEST_CASE("topology of the example meshes") {
  struct Expected {
    const char* name;
    long chi;
    std::size_t boundary_loops;
    long genus;
  };
  const Expected cases[] = {{"tetrahedron", 2, 0, 0}, {"cube", 2, 0, 0},  {"octahedron", 2, 0, 0},
                            {"saddle", 1, 1, 0},      {"torus", 0, 0, 1}, {"cylinder", 0, 2, 0},
                            {"l_prism", 2, 0, 0}};
  for (const auto& c : cases) {
    INFO(c.name);
    const auto t = topology(example(c.name));
    REQUIRE(t.euler_characteristic == c.chi);
    REQUIRE(t.boundary_loops == c.boundary_loops);
    REQUIRE(t.genus == c.genus);
    REQUIRE(t.components == 1);
    REQUIRE(t.is_disk() == (std::string(c.name) == "saddle"));
  }
}

TEST_CASE("Gauss-Bonnet holds on every example mesh") {
  for (const char* name : test::kExamples) {
    INFO(name);
    const auto mesh = example(name);
    REQUIRE(total_curvature(mesh) ==
            Approx(2 * pi * static_cast<double>(topology(mesh).euler_characteristic)).margin(1e-9));
  }
}

TEST_CASE("angle defects of the Platonic examples") {
  const auto tetrahedron = example("tetrahedron");  // a named mesh: the loop must not dangle
  for (const auto v : tetrahedron.vertices())
    REQUIRE(gaussian_curvature(tetrahedron, v) == Approx(pi));
  const auto cube = example("cube");
  for (const auto v : cube.vertices()) REQUIRE(gaussian_curvature(cube, v) == Approx(pi / 2));
  const auto octahedron = example("octahedron");
  for (const auto v : octahedron.vertices()) {
    REQUIRE(gaussian_curvature(octahedron, v) == Approx(2 * pi / 3));
  }
}

TEST_CASE("the saddle has negative curvature inside and boundary turning on its border") {
  const auto saddle = example("saddle");
  const auto center = SurfaceMesh::Vertex_index(24);  // (0, 0, 0) in the 7x7 grid
  REQUIRE(saddle.point(center) == origamizer::Point3(0, 0, 0));
  REQUIRE(gaussian_curvature(saddle, center) < 0);
  const auto corner = SurfaceMesh::Vertex_index(0);
  REQUIRE(saddle.is_border(corner));
}

TEST_CASE("boundary loops are vertex cycles along the border") {
  REQUIRE(origamizer::boundary_loops(example("cube")).empty());
  const auto saddle = origamizer::boundary_loops(example("saddle"));
  REQUIRE(saddle.size() == 1);
  REQUIRE(saddle[0].size() == 24);  // 6 x 4 boundary edges of the 7 x 7 grid
  const auto cylinder = origamizer::boundary_loops(example("cylinder"));
  REQUIRE(cylinder.size() == 2);
  for (const auto& loop : cylinder) REQUIRE(loop.size() == 12);
  // Consecutive loop vertices are joined by border edges.
  const auto mesh = example("saddle");
  for (std::size_t k = 0; k < saddle[0].size(); ++k) {
    const auto a = SurfaceMesh::Vertex_index(static_cast<SurfaceMesh::size_type>(saddle[0][k]));
    const auto b = SurfaceMesh::Vertex_index(
        static_cast<SurfaceMesh::size_type>(saddle[0][(k + 1) % saddle[0].size()]));
    const auto [h, found] = CGAL::halfedge(a, b, mesh);
    REQUIRE(found);
    REQUIRE(mesh.is_border(h));
  }
}

TEST_CASE("signed volumes and the clean side") {
  REQUIRE(origamizer::signed_volume(example("cube")) == Approx(1.0));
  REQUIRE(origamizer::signed_volume(example("tetrahedron")) == Approx(8.0 / 3.0));
  REQUIRE(origamizer::signed_volume(example("octahedron")) == Approx(4.0 / 3.0));

  auto inside_out = origamizer::read_obj(test::example_path("cube"));
  for (auto& face : inside_out.faces) std::reverse(face.begin(), face.end());
  const auto prepared = origamizer::prepare_mesh(inside_out);
  REQUIRE(prepared.flipped_to_clean_side);
  REQUIRE(origamizer::signed_volume(prepared.mesh) == Approx(1.0));

  auto saddle = example("saddle");
  REQUIRE_FALSE(origamizer::orient_clean_side_out(saddle));
}

TEST_CASE("simplify_triangle_mesh reduces faces and keeps the topology") {
  const auto saddle = example("saddle");
  const auto simplified = origamizer::simplify_triangle_mesh(saddle, 30);
  REQUIRE(simplified.number_of_faces() <= 30);
  REQUIRE(simplified.number_of_faces() > 0);
  REQUIRE(topology(simplified).is_disk());
  REQUIRE(origamizer::check_soup(origamizer::to_soup(simplified)).ok());
  REQUIRE_THROWS_WITH(origamizer::simplify_triangle_mesh(example("cube"), 4),
                      ContainsSubstring("needs a triangle mesh"));
}

TEST_CASE("prepare_mesh reports inconsistent orientation unless asked to fix it") {
  auto soup = origamizer::read_obj(test::example_path("cube"));
  std::reverse(soup.faces[2].begin(), soup.faces[2].end());
  REQUIRE_THROWS_WITH(origamizer::prepare_mesh(soup),
                      ContainsSubstring("not consistently oriented") &&
                          ContainsSubstring("make_orientation_consistent"));
  origamizer::PrepareOptions options;
  options.make_orientation_consistent = true;
  const auto prepared = origamizer::prepare_mesh(soup, options);
  REQUIRE(prepared.reoriented_faces == std::vector<std::size_t>{2});
  REQUIRE(origamizer::signed_volume(prepared.mesh) == Approx(1.0));
}
