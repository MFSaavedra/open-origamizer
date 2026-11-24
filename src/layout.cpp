// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

#include "origamizer/layout.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <numbers>

#include "detail/quadratic_program.hpp"
#include "origamizer/check.hpp"

namespace origamizer {
namespace {

using detail::kInfinity;
using detail::QuadraticProgram;
using std::numbers::pi;

std::size_t previous(std::size_t n, std::size_t size) { return (n + size - 1) % size; }

std::string format(double x) {
  char buffer[32];
  std::snprintf(buffer, sizeof buffer, "%.3g", x);
  return buffer;
}

}  // namespace

LinearEquations rotation_equations(const MoleculeMesh& mesh) {
  LinearEquations equations;
  equations.columns = mesh.molecules.size();
  for (std::size_t v = 0; v < mesh.stars.size(); ++v) {
    for (const StarSlot& slot : mesh.stars[v].slots) {
      equations.entries.push_back({v, slot.molecule, static_cast<double>(slot.sign)});
    }
    equations.rhs.push_back(mesh.stars[v].defect());
  }
  return equations;
}

std::vector<std::vector<double>> exterior_angles(const MoleculeMesh& mesh,
                                                 const std::vector<double>& theta) {
  ORIGAMIZER_CHECK(theta.size() == mesh.molecules.size(), "expected one angle per molecule");
  std::vector<std::vector<double>> exterior(mesh.stars.size());
  for (std::size_t v = 0; v < mesh.stars.size(); ++v) {
    const auto& slots = mesh.stars[v].slots;
    for (std::size_t n = 0; n < slots.size(); ++n) {
      const StarSlot& a = slots[previous(n, slots.size())];
      const StarSlot& b = slots[n];
      exterior[v].push_back(0.5 * a.sign * theta[a.molecule] + b.sector +
                            0.5 * b.sign * theta[b.molecule]);
    }
  }
  return exterior;
}

RotationResiduals rotation_residuals(const MoleculeMesh& mesh, const std::vector<double>& theta,
                                     double angle_margin) {
  ORIGAMIZER_CHECK(theta.size() == mesh.molecules.size(), "expected one angle per molecule");
  RotationResiduals residuals;
  const LinearEquations equations = rotation_equations(mesh);
  std::vector<double> lhs(equations.rhs.size(), 0);
  for (const SparseEntry& e : equations.entries) lhs[e.row] += e.value * theta[e.column];
  for (std::size_t k = 0; k < lhs.size(); ++k) {
    residuals.equations = std::max(residuals.equations, std::abs(lhs[k] - equations.rhs[k]));
  }
  for (std::size_t n = 0; n < mesh.molecules.size(); ++n) {
    const double violation = mesh.molecules[n].kind == Molecule::Kind::kBoundary
                                 ? pi - theta[n]                              // T4
                                 : std::abs(theta[n]) - (pi - angle_margin);  // T6
    residuals.inequalities = std::max(residuals.inequalities, violation);
  }
  for (const auto& star : exterior_angles(mesh, theta)) {
    for (const double angle : star) {  // T8
      residuals.inequalities =
          std::max({residuals.inequalities, -angle, angle - (pi - angle_margin)});
    }
  }
  return residuals;
}

RotationResult solve_rotation(const MoleculeMesh& mesh, const RotationOptions& options) {
  const double margin = options.angle_margin;
  ORIGAMIZER_CHECK(margin >= 0 && margin < pi / 4, "the angle margin must be in [0, pi/4)");
  RotationResult result;
  for (std::size_t v = 0; v < mesh.stars.size(); ++v) {
    const std::size_t slots = mesh.stars[v].slots.size();
    if (slots < 3) {
      result.message = "vertex " + std::to_string(v) + " has only " + std::to_string(slots) +
                       " molecules around it; T1 and T8 need at least three (enable the "
                       "boundary digons of TVCG10 Sec. 5.1)";
      return result;
    }
  }

  // Variables: one angle per molecule. Edge and digon molecules obey T6 and
  // enter the objective; boundary molecules obey T4.
  QuadraticProgram program;
  for (const Molecule& molecule : mesh.molecules) {
    const bool boundary = molecule.kind == Molecule::Kind::kBoundary;
    program.cost.push_back(0);
    program.diagonal.push_back(boundary ? 0 : 2);  // (1/2) * 2 * theta^2
    program.lower.push_back(boundary ? pi : -pi + margin);
    program.upper.push_back(boundary ? kInfinity : pi - margin);
  }
  const LinearEquations equations = rotation_equations(mesh);  // T15
  for (const double rhs : equations.rhs) program.add_row(rhs, rhs);
  for (const SparseEntry& e : equations.entries) {
    program.entries.push_back({e.row, e.column, e.value});
  }
  for (const VertexStar& star : mesh.stars) {  // T8: 0 <= Theta <= pi - margin
    for (std::size_t n = 0; n < star.slots.size(); ++n) {
      const StarSlot& a = star.slots[previous(n, star.slots.size())];
      const StarSlot& b = star.slots[n];
      const std::size_t row = program.add_row(-b.sector, pi - margin - b.sector);
      program.entries.push_back({row, a.molecule, 0.5 * a.sign});
      program.entries.push_back({row, b.molecule, 0.5 * b.sign});
    }
  }

  const detail::QpSolution solution = detail::solve(program);
  result.message = solution.message;
  if (solution.status != detail::SolveStatus::kOptimal) return result;
  const RotationResiduals residuals = rotation_residuals(mesh, solution.x, margin);
  ORIGAMIZER_CHECK(residuals.equations < 1e-6 && residuals.inequalities < 1e-6,
                   "the solver's angles violate the rotation constraints (equations by " +
                       format(residuals.equations) + ", inequalities by " +
                       format(residuals.inequalities) + ")");
  result.feasible = true;
  result.theta = solution.x;
  return result;
}

std::size_t numerical_rank(const LinearEquations& equations) {
  Eigen::MatrixXd a = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(equations.rhs.size()),
                                            static_cast<Eigen::Index>(equations.columns));
  for (const SparseEntry& e : equations.entries) {
    ORIGAMIZER_CHECK(e.row < equations.rhs.size() && e.column < equations.columns,
                     "a matrix entry is out of range");
    a(static_cast<Eigen::Index>(e.row), static_cast<Eigen::Index>(e.column)) += e.value;
  }
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(a);
  qr.setThreshold(1e-10);
  return static_cast<std::size_t>(qr.rank());
}

}  // namespace origamizer
