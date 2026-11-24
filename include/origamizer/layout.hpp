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

/// Mean length of the edge and digon molecules: the length scale of a layout,
/// used for default widths and relative tolerances.
double mean_molecule_length(const MoleculeMesh& mesh);

/// TVCG10 Eq. N6 (docs/tvcg2010-notes, Sec. 8.3), slot by slot: how much
/// wider the side of VTM(v) at slot n is than the molecule's narrower side,
/// delta(v, j_n) = 2 l max(0, -sin(theta(v, j_n) / 2)); 0 for a boundary
/// molecule, whose width variable is its side itself.
std::vector<std::vector<double>> width_excess(const MoleculeMesh& mesh,
                                              const std::vector<double>& theta);

/// T17: TVCG10 Eq. 2 for given angles, which is linear in the widths: two rows
/// per vertex (x, then y), one column per molecule (N6, N7). The side of
/// VTM(v) at slot n points along the sum of the exterior angles of slots 1 to n.
LinearEquations translation_equations(const MoleculeMesh& mesh, const std::vector<double>& theta);

/// The two objectives of the translation step (docs/notes.tex, M2).
enum class WidthObjective {
  /// A: minimize the sum of the widths, each at least the minimum width.
  kMinimumWidth,
  /// B: the paper's penalty. Phase 1 minimizes the total negative width; if
  /// it vanishes, phase 2 finds the nonnegative widths closest to a reference.
  kPaperPenalty,
};

struct TranslationOptions {
  WidthObjective objective = WidthObjective::kMinimumWidth;
  /// Objective A: the smallest allowed width. Negative selects the default,
  /// 1% of the mean molecule length.
  double minimum_width = -1;
  /// Objective B: the widths to stay close to, one per molecule. Empty selects
  /// the paper's starting point, the minimum-norm solution of T17.
  std::vector<double> reference;
};

struct TranslationResult {
  /// T17 holds and so do T5 and T7 (with the minimum width for objective A).
  bool feasible = false;
  /// One width per molecule: the narrower side of an edge or digon molecule,
  /// the boundary side w(i, o) of a boundary molecule. Objective B also
  /// returns phase 1's widths when T5 or T7 cannot hold; empty when T17 itself
  /// has no solution.
  std::vector<double> w;
  /// Objective B: phase 1's smallest total negative width (a length); 0 when
  /// T5 and T7 can hold.
  double violation = 0;
  std::string message;  ///< the solver's status, or why the step failed
};

/// The translation step (TVCG10 Sec. 6.3) for given angles, with either
/// objective (docs/notes.tex, M2). Checks the solution before returning it.
TranslationResult solve_translation(const MoleculeMesh& mesh, const std::vector<double>& theta,
                                    const TranslationOptions& options = {});

/// The largest residual of T17, and the smallest width (negative when T5 or
/// T7 fails).
struct TranslationResiduals {
  double equations = 0;
  double smallest_width = 0;
};

TranslationResiduals translation_residuals(const MoleculeMesh& mesh,
                                           const std::vector<double>& theta,
                                           const std::vector<double>& w);

/// The polygons of the molecule mesh placed in the plane (TVCG10 Sec. 3): a
/// congruent copy of every face of the disk and a segment for every digon.
/// The molecules fill the space between them.
struct Layout {
  /// For polygon p (MoleculeMesh numbering), the vertices of the disk at its
  /// corners: a face's vertices counterclockwise, or a digon's two ends.
  std::vector<std::vector<std::size_t>> vertices;
  /// The positions of those corners in the plane of the crease pattern.
  std::vector<std::vector<Point2>> corners;
  /// The largest distance between a side of a vertex-tucking molecule as
  /// placed and as TVCG10 Eq. 2 and 3 prescribe it.
  double closure_error = 0;
};

/// Places the polygons from the angles and widths (TVCG10 Sec. 3 and 4.1):
/// the first face at its own coordinates, then every polygon across an edge
/// or digon molecule from a placed one, by the geometry of the symmetric
/// edge-tucking molecule (docs/tvcg2010-notes, Fig. 1). Then checks every
/// side of every vertex-tucking molecule, which tests T15 and T17
/// independently of how they were assembled. Throws InvariantError if the
/// closure error exceeds 1e-6 of the mean molecule length.
Layout reconstruct_layout(const SurfaceMesh& disk, const MoleculeMesh& mesh,
                          const std::vector<double>& theta, const std::vector<double>& w);

/// Numerical rank of the matrix of a set of equations (column-pivoted QR, with
/// a threshold of 1e-10 relative to the largest pivot). Used for the
/// degree-of-freedom check of TVCG10 Sec. 6.1 (docs/tvcg2010-notes, Open
/// question 5).
std::size_t numerical_rank(const LinearEquations& equations);

}  // namespace origamizer
