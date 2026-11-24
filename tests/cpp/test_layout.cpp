// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

#include <Eigen/Dense>
#include <array>
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

origamizer::SurfaceMesh disk_of(const std::string& name) {
  const auto mesh = origamizer::prepare_mesh(origamizer::read_obj(test::example_path(name))).mesh;
  return origamizer::cut_to_disk(mesh).disk;
}

origamizer::MoleculeMesh molecules_of(const std::string& name,
                                      const origamizer::MoleculeOptions& options = {}) {
  return origamizer::build_molecule_mesh(disk_of(name), options);
}

/// Examples whose rotation-step angles admit nonnegative widths; on the cube
/// and the L prism they do not (TVCG10 Sec. 6.4), which needs remapping.
const std::array<const char*, 5> kTranslatable = {"tetrahedron", "octahedron", "saddle", "torus",
                                                  "cylinder"};

Eigen::MatrixXd dense(const origamizer::LinearEquations& equations) {
  Eigen::MatrixXd a = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(equations.rhs.size()),
                                            static_cast<Eigen::Index>(equations.columns));
  for (const auto& e : equations.entries) {
    a(static_cast<Eigen::Index>(e.row), static_cast<Eigen::Index>(e.column)) += e.value;
  }
  return a;
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

TEST_CASE("an infeasible rotation step reports its least violation") {
  // One synthetic star, sectors alternating 3.0 and 0.05: T8 at the narrow
  // corners needs theta_{m-1} + theta_m >= -0.1 for three disjoint pairs, so
  // the sum of theta is at least -0.3, while T1 asks for 2 pi - 9.15. The
  // least total violation is -(2 pi - 9.15) / 2 - 3 * 0.05.
  origamizer::MoleculeMesh mesh;
  mesh.faces = 6;
  origamizer::VertexStar star;
  for (std::size_t n = 0; n < 6; ++n) {
    mesh.molecules.push_back({Molecule::Kind::kEdge, 0, n + 1, 1.0});
    star.slots.push_back({n, 1, n % 2 == 0 ? 3.0 : 0.05, n});
  }
  mesh.stars.push_back(star);
  mesh.edge_molecules = 6;
  const auto rotation = origamizer::solve_rotation(mesh);
  REQUIRE_FALSE(rotation.feasible);
  REQUIRE(rotation.theta.size() == 6);
  const double expected = -(2 * pi - 9.15) / 2 - 0.15;
  REQUIRE(rotation.violation == Approx(expected).epsilon(1e-9));
  REQUIRE_THAT(rotation.message, ContainsSubstring("least total violation"));
  double sum = 0;
  for (const auto& violation : rotation.violations) {
    REQUIRE(violation.kind == origamizer::RotationViolation::Kind::kT8);
    REQUIRE(violation.slot % 2 == 1);  // the narrow corners
    sum += violation.amount;
  }
  REQUIRE(sum == Approx(rotation.violation));
  // T15 still holds for the least-violation angles.
  double total = 0;
  for (const double t : rotation.theta) total += t;
  REQUIRE(total == Approx(2 * pi - 9.15).margin(1e-12));
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

TEST_CASE("translation step, objective A, where nonnegative widths exist") {
  for (const char* name : kTranslatable) {
    INFO(name);
    const auto mesh = molecules_of(name);
    const auto rotation = origamizer::solve_rotation(mesh);
    REQUIRE(rotation.feasible);
    const auto translation = origamizer::solve_translation(mesh, rotation.theta);
    INFO(translation.message);
    REQUIRE(translation.feasible);
    const double scale = origamizer::mean_molecule_length(mesh);
    const auto r = origamizer::translation_residuals(mesh, rotation.theta, translation.w);
    REQUIRE(r.equations < 1e-9 * scale);
    REQUIRE(r.smallest_width > 0.01 * scale - 1e-9 * scale);
  }
}

TEST_CASE("objective A honors a requested minimum width") {
  const auto mesh = molecules_of("octahedron");
  const auto theta = origamizer::solve_rotation(mesh).theta;
  origamizer::TranslationOptions options;
  options.minimum_width = 0.1;
  const auto translation = origamizer::solve_translation(mesh, theta, options);
  REQUIRE(translation.feasible);
  REQUIRE(origamizer::translation_residuals(mesh, theta, translation.w).smallest_width >
          0.1 - 1e-9);
}

TEST_CASE("objective B finds the nonnegative widths closest to the minimum-norm solution") {
  for (const char* name : kTranslatable) {
    INFO(name);
    const auto mesh = molecules_of(name);
    const auto theta = origamizer::solve_rotation(mesh).theta;
    origamizer::TranslationOptions options;
    options.objective = origamizer::WidthObjective::kPaperPenalty;
    const auto b = origamizer::solve_translation(mesh, theta, options);
    INFO(b.message);
    REQUIRE(b.feasible);
    REQUIRE(b.violation == 0);
    const double scale = origamizer::mean_molecule_length(mesh);
    const auto r = origamizer::translation_residuals(mesh, theta, b.w);
    REQUIRE(r.equations < 1e-9 * scale);
    REQUIRE(r.smallest_width > -1e-12 * scale);
    // Objective A's widths are also nonnegative solutions of T17, so they are
    // at least as far from the reference C_w^+ b.
    const auto equations = origamizer::translation_equations(mesh, theta);
    const Eigen::VectorXd rhs = Eigen::Map<const Eigen::VectorXd>(
        equations.rhs.data(), static_cast<Eigen::Index>(equations.rhs.size()));
    const Eigen::VectorXd reference = dense(equations).completeOrthogonalDecomposition().solve(rhs);
    const auto a = origamizer::solve_translation(mesh, theta);
    const Eigen::VectorXd wa = Eigen::Map<const Eigen::VectorXd>(a.w.data(), reference.size());
    const Eigen::VectorXd wb = Eigen::Map<const Eigen::VectorXd>(b.w.data(), reference.size());
    REQUIRE((wb - reference).norm() <= (wa - reference).norm() + 1e-9 * scale);
  }
}

TEST_CASE("objective B stays at a feasible reference") {
  const auto mesh = molecules_of("torus");
  const auto theta = origamizer::solve_rotation(mesh).theta;
  origamizer::TranslationOptions options;
  options.objective = origamizer::WidthObjective::kPaperPenalty;
  options.reference = origamizer::solve_translation(mesh, theta).w;  // objective A's widths
  const auto b = origamizer::solve_translation(mesh, theta, options);
  REQUIRE(b.feasible);
  for (std::size_t n = 0; n < b.w.size(); ++n) {
    REQUIRE(b.w[n] == Approx(options.reference[n]).margin(1e-9));
  }
}

TEST_CASE("the cube and the L prism need negative widths with the rotation step's angles") {
  // TVCG10 Sec. 6.4: fixing theta first can leave T17 without a nonnegative
  // solution. Before the solutions were polished, solver noise of 1e-8 in
  // theta let objective A report the cube feasible with widths near 1e8.
  for (const char* name : {"cube", "l_prism"}) {
    INFO(name);
    const auto mesh = molecules_of(name);
    const auto theta = origamizer::solve_rotation(mesh).theta;
    const auto a = origamizer::solve_translation(mesh, theta);
    REQUIRE_FALSE(a.feasible);
    REQUIRE_THAT(a.message, ContainsSubstring("objective A: Infeasible"));
    REQUIRE(a.w.empty());
    origamizer::TranslationOptions options;
    options.objective = origamizer::WidthObjective::kPaperPenalty;
    const auto b = origamizer::solve_translation(mesh, theta, options);
    REQUIRE_FALSE(b.feasible);
    REQUIRE(b.violation > 1.0);
    REQUIRE_THAT(b.message, ContainsSubstring("T5 and T7 cannot hold"));
    const double scale = origamizer::mean_molecule_length(mesh);
    const auto r = origamizer::translation_residuals(mesh, theta, b.w);
    REQUIRE(r.equations < 1e-9 * scale);
    REQUIRE(r.smallest_width < 0);
    double negative = 0;
    for (const double w : b.w) negative += std::max(0.0, -w);
    REQUIRE(negative == Approx(b.violation).epsilon(1e-9));
  }
}

TEST_CASE("T17 is dependent but consistent on four examples") {
  // Rank 2V - 1 on the tetrahedron, octahedron, torus, and cylinder (with b in
  // the range, since the translation step succeeds there); full row rank on
  // the others. The cylinder's dependency is exact only with polished angles:
  // with HiGHS's own, its smallest singular value was 3e-9.
  for (const char* name : test::kExamples) {
    INFO(name);
    const auto mesh = molecules_of(name);
    const auto theta = origamizer::solve_rotation(mesh).theta;
    const std::string n = name;
    const bool dependent =
        n == "tetrahedron" || n == "octahedron" || n == "torus" || n == "cylinder";
    REQUIRE(origamizer::numerical_rank(origamizer::translation_equations(mesh, theta)) ==
            2 * mesh.stars.size() - (dependent ? 1 : 0));
  }
}

TEST_CASE("the layout closes for widths that solve T17, negative or not") {
  // Objective B returns phase 1's widths when nonnegative ones do not exist;
  // they still solve T17, so the molecule mesh must close either way.
  for (const char* name : test::kExamples) {
    INFO(name);
    const auto disk = disk_of(name);
    const auto mesh = origamizer::build_molecule_mesh(disk);
    const auto theta = origamizer::solve_rotation(mesh).theta;
    origamizer::TranslationOptions options;
    options.objective = origamizer::WidthObjective::kPaperPenalty;
    const auto translation = origamizer::solve_translation(mesh, theta, options);
    INFO(translation.message);
    REQUIRE(translation.w.size() == mesh.molecules.size());
    const auto layout = origamizer::reconstruct_layout(disk, mesh, theta, translation.w);
    REQUIRE(layout.closure_error < 1e-9 * origamizer::mean_molecule_length(mesh));
    REQUIRE(layout.corners.size() == mesh.polygons());
    // Faces are placed congruently: edge lengths as on the surface.
    for (const auto f : disk.faces()) {
      const auto p = static_cast<std::size_t>(f);
      const auto& corners = layout.corners[p];
      std::size_t k = 0;
      for (const auto h : CGAL::halfedges_around_face(disk.halfedge(f), disk)) {
        const double placed =
            std::sqrt(CGAL::squared_distance(corners[k], corners[(k + 1) % corners.size()]));
        const double original = std::sqrt(CGAL::squared_distance(
            disk.point(disk.source(disk.next(h))), disk.point(disk.target(disk.next(h)))));
        REQUIRE(placed == Approx(original).margin(1e-12));
        ++k;
      }
    }
  }
}

TEST_CASE("the reconstruction rejects widths that do not solve T17") {
  const auto disk = disk_of("octahedron");
  const auto mesh = origamizer::build_molecule_mesh(disk);
  const auto theta = origamizer::solve_rotation(mesh).theta;
  auto w = origamizer::solve_translation(mesh, theta).w;
  w[0] += 0.1;
  REQUIRE_THROWS_WITH(origamizer::reconstruct_layout(disk, mesh, theta, w),
                      ContainsSubstring("does not close"));
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
