// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <numbers>
#include <string>

#include "origamizer/cut.hpp"
#include "origamizer/layout.hpp"
#include "support.hpp"

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;
using origamizer::Molecule;
using std::numbers::pi;

namespace {

origamizer::MoleculeMesh molecules_of(const std::string& name,
                                      const origamizer::MoleculeOptions& options = {}) {
  const auto mesh = origamizer::prepare_mesh(origamizer::read_obj(test::example_path(name))).mesh;
  return origamizer::build_molecule_mesh(origamizer::cut_to_disk(mesh).disk, options);
}

}  // namespace

TEST_CASE("the rotation step is feasible on every example with boundary digons") {
  for (const char* name : test::kExamples) {
    INFO(name);
    const auto mesh = molecules_of(name);
    const auto rotation = origamizer::solve_rotation(mesh);
    INFO(rotation.message);
    REQUIRE(rotation.feasible);
    const auto residuals = origamizer::rotation_residuals(mesh, rotation.theta, 1e-3);
    REQUIRE(residuals.equations < 1e-6);
    REQUIRE(residuals.inequalities < 1e-6);
    // Gauss-Bonnet (docs/tvcg2010-notes, Sec. 4): interior molecules cancel in
    // the sum of T1, so the boundary turning adds up to one full turn.
    double turning = 0;
    for (std::size_t n = 0; n < mesh.molecules.size(); ++n) {
      if (mesh.molecules[n].kind == Molecule::Kind::kBoundary) turning += rotation.theta[n] - pi;
    }
    REQUIRE(turning == Approx(2 * pi).margin(1e-6));
  }
}

TEST_CASE("without boundary digons, a vertex with fewer than three molecules is reported") {
  const auto mesh = molecules_of("cube", {.boundary_digons = false});
  const auto rotation = origamizer::solve_rotation(mesh);
  REQUIRE_FALSE(rotation.feasible);
  REQUIRE(rotation.theta.empty());
  REQUIRE_THAT(rotation.message, ContainsSubstring("need at least three"));
}

TEST_CASE("exterior angles of a star sum to 2 pi when T1 holds") {
  const auto mesh = molecules_of("octahedron");
  const auto rotation = origamizer::solve_rotation(mesh);
  REQUIRE(rotation.feasible);
  for (const auto& star : origamizer::exterior_angles(mesh, rotation.theta)) {
    double sum = 0;
    for (const double angle : star) {
      REQUIRE(angle >= -1e-9);
      REQUIRE(angle <= pi - 1e-3 + 1e-9);
      sum += angle;
    }
    REQUIRE(sum == Approx(2 * pi));
  }
}

TEST_CASE("the rotation equations have full row rank on a disk") {
  // One row per vertex: an incidence matrix of a connected graph (rank V - 1)
  // plus the boundary molecules' columns, which add the missing rank.
  for (const char* name : {"tetrahedron", "cube", "saddle", "torus"}) {
    INFO(name);
    const auto mesh = molecules_of(name);
    REQUIRE(origamizer::numerical_rank(origamizer::rotation_equations(mesh)) == mesh.stars.size());
  }
}
