// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

/// \file
/// The two-step mapping of TVCG10 Sec. 6 (milestone M2). The rotation step
/// chooses the molecule angles theta; the translation step chooses the widths
/// w. Both are linear or quadratic programs solved with HiGHS (docs/notes.tex,
/// M2). Equation tags refer to docs/tvcg2010-notes.

#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "origamizer/molecules.hpp"

namespace origamizer {

/// A nonzero of a sparse matrix.
struct SparseEntry {
  std::size_t row = 0;
  std::size_t column = 0;
  double value = 0;
};

/// Linear equations A x = rhs with A sparse; repeated entries add up.
struct LinearEquations {
  std::size_t columns = 0;
  std::vector<SparseEntry> entries;
  std::vector<double> rhs;
};

/// T15: TVCG10 Eq. 1 at every vertex, one row per vertex and one column per
/// molecule: the sum of theta(v, j_n) over the star equals its defect.
LinearEquations rotation_equations(const MoleculeMesh& mesh);

/// TVCG10 Eq. 3: the exterior angle of every vertex-tucking molecule at each
/// slot, Theta[v][n] = theta(v, j_{n-1}) / 2 + alpha(v, j_n) + theta(v, j_n) / 2.
std::vector<std::vector<double>> exterior_angles(const MoleculeMesh& mesh,
                                                 const std::vector<double>& theta);

struct RotationOptions {
  /// Margin of the strict inequalities T6 and T8, in radians (docs/notes.tex, M2).
  double angle_margin = 1e-3;
};

struct RotationResult {
  bool feasible = false;
  /// One angle per molecule: theta(i, j) for its orientation i -> j, or
  /// theta(i, o) for a boundary molecule. Empty unless feasible.
  std::vector<double> theta;
  std::string message;  ///< the solver's status, or why the step cannot succeed
};

/// The largest violation of T15, and of T4, T6, and T8 with the margin
/// (0 when they all hold).
struct RotationResiduals {
  double equations = 0;
  double inequalities = 0;
};

RotationResiduals rotation_residuals(const MoleculeMesh& mesh, const std::vector<double>& theta,
                                     double angle_margin);

/// The rotation step (TVCG10 Sec. 6.2) as a quadratic program: minimize the
/// sum of theta^2 over the edge and digon molecules subject to T15, T4, T6,
/// and T8, with T6 and T8 tightened by the margin (docs/notes.tex, M2). A
/// star with fewer than three slots makes T1 and T8 contradict each other;
/// the result then says so without calling the solver. The solution is
/// checked against the constraints before it is returned.
RotationResult solve_rotation(const MoleculeMesh& mesh, const RotationOptions& options = {});

/// Numerical rank of the matrix of a set of equations (column-pivoted QR, with
/// a threshold of 1e-10 relative to the largest pivot). Used for the
/// degree-of-freedom check of TVCG10 Sec. 6.1 (docs/tvcg2010-notes, Open
/// question 5).
std::size_t numerical_rank(const LinearEquations& equations);

}  // namespace origamizer
