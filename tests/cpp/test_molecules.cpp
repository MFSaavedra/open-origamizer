// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <numbers>
#include <string>

#include "origamizer/check.hpp"
#include "origamizer/cut.hpp"
#include "origamizer/molecules.hpp"
#include "support.hpp"

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;
using origamizer::build_molecule_mesh;
using origamizer::Molecule;
using std::numbers::pi;

namespace {

origamizer::SurfaceMesh disk_of(const origamizer::PolygonSoup& soup) {
  return origamizer::cut_to_disk(origamizer::prepare_mesh(soup).mesh).disk;
}

/// The vertex at the other end of slot `slot` of the star of v.
std::size_t neighbor(const origamizer::MoleculeMesh& mesh, std::size_t v, std::size_t slot) {
  const Molecule& m = mesh.molecules[mesh.stars[v].slots[slot].molecule];
  return m.i == v ? m.j : m.i;
}

}  // namespace

namespace {
const origamizer::MoleculeOptions kNoDigons{.boundary_digons = false};
}

TEST_CASE("a single triangle: one slot per corner without digons, three with them") {
  const auto disk = disk_of(test::flat_polygon({{0, 0}, {1, 0}, {0.5, std::sqrt(3.0) / 2}}));
  const auto bare = build_molecule_mesh(disk, kNoDigons);
  REQUIRE(bare.edge_molecules == 0);
  REQUIRE(bare.digon_molecules == 0);
  REQUIRE(bare.boundary_molecules() == 3);
  for (const auto& star : bare.stars) {
    REQUIRE(star.boundary);
    REQUIRE(star.slots.size() == 1);
    REQUIRE(star.slots[0].sector == Approx(pi / 3));
    REQUIRE(star.defect() == Approx(2 * pi - pi / 3));
  }

  // TVCG10 Sec. 5.1: a digon molecule per boundary edge; each corner star
  // becomes [digon (corner pi/3), outside (0), digon (0)].
  const auto mesh = build_molecule_mesh(disk);
  REQUIRE(mesh.digon_molecules == 3);
  REQUIRE(mesh.boundary_molecules() == 3);
  for (const auto& star : mesh.stars) {
    REQUIRE(star.slots.size() == 3);
    REQUIRE(mesh.molecules[star.slots[0].molecule].kind == Molecule::Kind::kDigon);
    REQUIRE(mesh.molecules[star.slots[1].molecule].kind == Molecule::Kind::kBoundary);
    REQUIRE(mesh.molecules[star.slots[2].molecule].kind == Molecule::Kind::kDigon);
    REQUIRE(star.slots[0].sector == Approx(pi / 3));
    REQUIRE(star.slots[1].sector == 0);
    REQUIRE(star.slots[2].sector == 0);
    REQUIRE(star.defect() == Approx(2 * pi - pi / 3));
  }
}

TEST_CASE("a split square has one edge molecule shared by two stars with opposite signs") {
  // Square 0-1-2-3 split along the diagonal 0-2.
  const auto soup =
      test::make_soup({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, {{0, 1, 2}, {0, 2, 3}});
  const auto with_digons = build_molecule_mesh(disk_of(soup));
  REQUIRE(with_digons.digon_molecules == 4);
  REQUIRE(with_digons.stars[0].slots.size() == 4);
  REQUIRE(with_digons.stars[1].slots.size() == 3);
  const auto mesh = build_molecule_mesh(disk_of(soup), kNoDigons);
  REQUIRE(mesh.edge_molecules == 1);
  REQUIRE(mesh.boundary_molecules() == 4);
  REQUIRE(mesh.molecules[0].i == 0);
  REQUIRE(mesh.molecules[0].j == 2);
  REQUIRE(mesh.molecules[0].length == Approx(std::sqrt(2.0)));
  REQUIRE(mesh.stars[0].slots.size() == 2);
  REQUIRE(mesh.stars[2].slots.size() == 2);
  REQUIRE(mesh.stars[1].slots.size() == 1);
  int sign_sum = 0;
  for (const std::size_t v : {std::size_t{0}, std::size_t{2}}) {
    double sectors = 0;
    for (const auto& slot : mesh.stars[v].slots) {
      sectors += slot.sector;
      if (slot.molecule == 0) sign_sum += slot.sign;
    }
    REQUIRE(sectors == Approx(pi / 2));
  }
  REQUIRE(sign_sum == 0);
}

TEST_CASE("star slots run counterclockwise around an interior vertex") {
  // A flat hexagonal fan around vertex 0, faces counterclockwise.
  origamizer::PolygonSoup soup;
  soup.points.emplace_back(0, 0, 0);
  for (int k = 0; k < 6; ++k) {
    soup.points.emplace_back(std::cos(k * pi / 3), std::sin(k * pi / 3), 0);
  }
  for (std::size_t k = 0; k < 6; ++k) soup.faces.push_back({0, 1 + k, 1 + (k + 1) % 6});
  const auto disk = disk_of(soup);
  const auto mesh = build_molecule_mesh(disk);
  const auto& star = mesh.stars[0];
  REQUIRE_FALSE(star.boundary);
  REQUIRE(star.slots.size() == 6);
  // Each next slot is 60 degrees further counterclockwise.
  const auto angle_of = [&](std::size_t slot) {
    const auto p = disk.point(origamizer::SurfaceMesh::Vertex_index(
        static_cast<origamizer::SurfaceMesh::size_type>(neighbor(mesh, 0, slot))));
    return std::atan2(p.y(), p.x());
  };
  for (std::size_t n = 0; n < 6; ++n) {
    double step = angle_of((n + 1) % 6) - angle_of(n);
    if (step < 0) step += 2 * pi;
    REQUIRE(step == Approx(pi / 3));
    REQUIRE(star.slots[n].sector == Approx(pi / 3));
  }
  REQUIRE(star.defect() == Approx(0.0).margin(1e-12));
}

TEST_CASE("defects of every example disk satisfy Gauss-Bonnet") {
  // Summing T1's right side over a disk gives 2 pi plus pi per boundary vertex
  // (docs/tvcg2010-notes, Sec. 4, consistency check).
  for (const char* name : test::kExamples) {
    INFO(name);
    const auto disk = disk_of(origamizer::read_obj(test::example_path(name)));
    const auto mesh = build_molecule_mesh(disk);
    double total = 0;
    std::size_t boundary_vertices = 0;
    for (const auto& star : mesh.stars) {
      total += star.defect();
      boundary_vertices += star.boundary;
    }
    REQUIRE(mesh.boundary_molecules() == boundary_vertices);
    REQUIRE(total == Approx(2 * pi + pi * static_cast<double>(boundary_vertices)).margin(1e-9));
    std::size_t interior_edges = 0, boundary_edges = 0;
    for (const auto e : disk.edges()) (disk.is_border(e) ? boundary_edges : interior_edges) += 1;
    REQUIRE(mesh.edge_molecules == interior_edges);
    REQUIRE(mesh.digon_molecules == boundary_edges);
    // T8 needs at least three slots around every vertex; the digons provide them.
    for (const auto& star : mesh.stars) REQUIRE(star.slots.size() >= 3);
  }
}

TEST_CASE("build_molecule_mesh rejects a surface that is not a disk") {
  const auto cube = origamizer::prepare_mesh(origamizer::read_obj(test::example_path("cube"))).mesh;
  REQUIRE_THROWS_WITH(build_molecule_mesh(cube), ContainsSubstring("cut the surface first"));
}
