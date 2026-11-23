// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <map>
#include <numbers>
#include <string>

#include "origamizer/check.hpp"
#include "origamizer/cut.hpp"
#include "origamizer/mesh.hpp"
#include "support.hpp"

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;
using origamizer::cut_along;
using origamizer::cut_to_disk;
using origamizer::IndexPair;
using origamizer::SurfaceMesh;
using origamizer::topology;

namespace {

SurfaceMesh example(const std::string& name) {
  return origamizer::prepare_mesh(origamizer::read_obj(test::example_path(name))).mesh;
}

/// Degree of each vertex in the graph formed by the cut edges.
std::map<std::size_t, int> cut_degrees(const std::vector<IndexPair>& edges) {
  std::map<std::size_t, int> degree;
  for (const auto& [a, b] : edges) {
    ++degree[a];
    ++degree[b];
  }
  return degree;
}

/// Invariants every cut must satisfy, whatever edges were chosen.
void require_valid_cut(const SurfaceMesh& mesh, const origamizer::CutResult& cut) {
  const auto t = topology(cut.disk);
  REQUIRE(t.is_disk());
  REQUIRE(cut.disk.number_of_faces() == mesh.number_of_faces());
  REQUIRE(cut.original_vertex.size() == cut.disk.number_of_vertices());
  const auto disk_soup = origamizer::to_soup(cut.disk);
  const auto mesh_soup = origamizer::to_soup(mesh);
  for (std::size_t v = 0; v < disk_soup.points.size(); ++v) {
    REQUIRE(disk_soup.points[v] == mesh_soup.points[cut.original_vertex[v]]);
  }
  for (std::size_t f = 0; f < disk_soup.faces.size(); ++f) {
    std::vector<std::size_t> mapped;
    for (const std::size_t v : disk_soup.faces[f]) mapped.push_back(cut.original_vertex[v]);
    const auto& original = mesh_soup.faces[f];
    REQUIRE(mapped.size() == original.size());
    const auto start = std::find(original.begin(), original.end(), mapped[0]);
    REQUIRE(start != original.end());
    for (std::size_t k = 0; k < mapped.size(); ++k) {
      REQUIRE(mapped[k] ==
              original[(static_cast<std::size_t>(start - original.begin()) + k) % original.size()]);
    }
  }
  // Each cut edge adds two boundary edges to the original boundary.
  REQUIRE(test::boundary_length(cut.disk) ==
          test::boundary_length(mesh) + 2 * cut.cut_edges.size());
  // A disk has total curvature 2 pi (Gauss-Bonnet with chi = 1).
  double total = 0;
  for (const auto v : cut.disk.vertices()) total += origamizer::gaussian_curvature(cut.disk, v);
  REQUIRE(total == Approx(2 * std::numbers::pi).margin(1e-9));
}

}  // namespace

TEST_CASE("every example mesh cuts into a disk") {
  for (const char* name : test::kExamples) {
    INFO(name);
    const auto mesh = example(name);
    require_valid_cut(mesh, cut_to_disk(mesh));
  }
}

TEST_CASE("closed genus-0 meshes are cut along one path between their most curved vertices") {
  // All tetrahedron vertices are adjacent, so the one-edge path 0-1 is extended to 0-1-2.
  const auto tetrahedron = cut_to_disk(example("tetrahedron"));
  REQUIRE(tetrahedron.cut_edges == std::vector<IndexPair>{{0, 1}, {1, 2}});
  REQUIRE(tetrahedron.disk.number_of_vertices() == 5);  // vertex 1 split; the ends stay whole

  const auto cube = cut_to_disk(example("cube"));
  REQUIRE(cube.cut_edges.size() == 3);  // from vertex 0 to the opposite corner, vertex 7
  const auto degree = cut_degrees(cube.cut_edges);
  REQUIRE(degree.at(0) == 1);
  REQUIRE(degree.at(7) == 1);
  REQUIRE(cube.disk.number_of_vertices() == 10);

  const auto octahedron = cut_to_disk(example("octahedron"));
  REQUIRE(octahedron.cut_edges.size() == 2);  // between opposite vertices 0 and 1
  REQUIRE(cut_degrees(octahedron.cut_edges).at(1) == 1);
}

TEST_CASE("a disk needs no cut") {
  const auto saddle = example("saddle");
  const auto cut = cut_to_disk(saddle);
  REQUIRE(cut.cut_edges.empty());
  REQUIRE(cut.disk.number_of_vertices() == saddle.number_of_vertices());
}

TEST_CASE("a torus is cut along loops, with no dangling cut edges") {
  const auto torus = example("torus");
  const auto cut = cut_to_disk(torus);
  REQUIRE(cut.cut_edges.size() >= 2);
  for (const auto& [v, d] : cut_degrees(cut.cut_edges)) {
    INFO("vertex " << v);
    REQUIRE(d >= 2);
  }
  // A different root gives a different, equally valid cut.
  origamizer::CutOptions options;
  options.root_face = 50;
  require_valid_cut(torus, cut_to_disk(torus, options));
}

TEST_CASE("an open cylinder is cut along one path joining its two boundary loops") {
  const auto cylinder = example("cylinder");
  const auto cut = cut_to_disk(cylinder);
  REQUIRE(cut.cut_edges.size() == 1);  // one rung between the rings
  const auto [a, b] = cut.cut_edges[0];
  REQUIRE((a < 12) != (b < 12));  // one end on each boundary loop
}

TEST_CASE("explicit cut edges are checked") {
  const auto tetrahedron = example("tetrahedron");
  const auto path = cut_along(tetrahedron, {{0, 1}, {1, 2}});
  REQUIRE(path.disk.number_of_vertices() == 5);  // vertex 1 is split in two
  REQUIRE_THROWS_WITH(cut_along(tetrahedron, {{0, 1}, {1, 2}, {2, 0}}),
                      ContainsSubstring("does not leave one surface"));
  REQUIRE_THROWS_WITH(cut_along(tetrahedron, {{0, 7}}),
                      ContainsSubstring("does not name two vertices"));
  REQUIRE_THROWS_WITH(cut_along(tetrahedron, {{2, 3}}), ContainsSubstring("would be glued back"));

  const auto octahedron = example("octahedron");
  // The equator 2-0-3-1 separates the two pyramids.
  REQUIRE_THROWS_WITH(cut_along(octahedron, {{0, 2}, {0, 3}, {1, 2}, {1, 3}}),
                      ContainsSubstring("does not leave one surface"));
  REQUIRE_THROWS_WITH(cut_along(octahedron, {{0, 1}}), ContainsSubstring("not an edge"));

  const auto saddle = example("saddle");
  REQUIRE_THROWS_WITH(cut_along(saddle, {{0, 1}}), ContainsSubstring("boundary edge"));
  // Cutting the torus along one loop only leaves an annulus.
  const auto torus = example("torus");
  std::vector<IndexPair> one_loop;
  for (std::size_t i = 0; i < 12; ++i) {
    const std::size_t a = i * 8, b = ((i + 1) % 12) * 8;
    one_loop.push_back({std::min(a, b), std::max(a, b)});
  }
  REQUIRE_THROWS_WITH(cut_along(torus, one_loop), ContainsSubstring("not a topological disk"));
}

TEST_CASE("cut_to_disk rejects a root face that does not exist") {
  origamizer::CutOptions options;
  options.root_face = 1000;
  REQUIRE_THROWS_WITH(cut_to_disk(example("cube"), options), ContainsSubstring("root face 1000"));
}
