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

#include "origamizer/cut.hpp"
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
  /// Angles to stay close to, one per molecule: the objective becomes the
  /// sum of (theta - reference)^2 over the edge and digon molecules. Empty
  /// means zero, the least-norm angles. Remapping (remap_angles) supplies one.
  std::vector<double> reference;
};

/// An inequality of the rotation step that the least-violation angles break.
struct RotationViolation {
  enum class Kind {
    kT4,  ///< theta(i, o) >= pi at a boundary molecule
    kT6,  ///< |theta| <= pi - margin at an edge or digon molecule
    kT8,  ///< 0 <= Theta <= pi - margin at a corner of a vertex-tucking molecule
  };
  Kind kind = Kind::kT6;
  std::size_t molecule = 0;  ///< for T4 and T6
  std::size_t vertex = 0;    ///< for T8: the vertex, and the slot whose corner it is
  std::size_t slot = 0;
  double amount = 0;  ///< by how much, in radians
};

struct RotationResult {
  bool feasible = false;
  /// One angle per molecule: theta(i, j) for its orientation i -> j, or
  /// theta(i, o) for a boundary molecule. When the step is infeasible, the
  /// angles of least total violation (T15 still holds), or empty if a star
  /// has fewer than three slots.
  std::vector<double> theta;
  /// When infeasible: the least total violation of T4, T6, and T8 (radians),
  /// and the inequalities those angles break, largest first.
  double violation = 0;
  std::vector<RotationViolation> violations;
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
/// the result then says so without calling the solver. If the program is
/// infeasible, an elastic linear program keeps T15 and minimizes the total
/// violation of the inequalities, to report which ones fail. The solution is
/// checked against the constraints before it is returned.
RotationResult solve_rotation(const MoleculeMesh& mesh, const RotationOptions& options = {});

/// Remapping (TVCG10 Sec. 6.5): reference angles for the rotation step from a
/// crossing-free drawing of the molecule graph, whose nodes are the polygons
/// (faces and digons), whose edges are the molecules, and whose faces are the
/// vertex-tucking molecules. Draws the graph by Tutte's barycentric method
/// with the digons evenly on a circle, measures the exterior angle at every
/// corner of every vertex-tucking molecule, and fits the angles to them by
/// least squares through T3. Needs the boundary digons (docs/notes.tex, M2).
std::vector<double> remap_angles(const MoleculeMesh& mesh);

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
  /// The boundary of the paper, counterclockwise: for each boundary vertex in
  /// order along the boundary, its copies in the digon before and after its
  /// boundary molecule. Between those two runs the boundary molecule's side;
  /// from one boundary vertex to the next, a digon. Empty without digons.
  std::vector<Point2> outline;
  /// Every edge-tucking molecule as a quadrilateral A0 B0 B1 A1 (the legs in
  /// the polygons before and after it, counterclockwise around its vertex
  /// i), indexed by molecule; empty for boundary molecules.
  std::vector<std::vector<Point2>> edge_molecules;
  /// Every vertex-tucking molecule, its corners counterclockwise.
  std::vector<std::vector<Point2>> vertex_molecules;
};

/// Whether a layout is a crease pattern on convex paper without overlap
/// (TVCG10 Sec. 4.2.1): with a convex outline, the molecule mesh cannot
/// overlap unless a cell intersects itself or flips, so it is enough that
/// every molecule is a convex polygon with positive area or degenerate (zero
/// area, as zero widths allow: a point or a segment, which covers nothing).
/// T4 to T8 guarantee this; the check verifies it on the placed polygons.
struct OverlapReport {
  /// The smallest signed area of an edge- or vertex-tucking molecule.
  double smallest_area = 0;
  /// Molecules of zero area (within 1e-9 of the squared mean edge length).
  std::size_t degenerate_molecules = 0;
  /// The most negative turning angle at a corner of a molecule of positive
  /// area or of the outline (a reflex corner), or 0 if there is none.
  /// Corners on sides of zero length are skipped.
  double reflex_turn = 0;
  /// The total turning of the outline: 2 pi for a simple convex outline.
  double outline_turning = 0;
  /// The areas of the faces and molecules minus the area of the outline,
  /// relative to the outline's area (Green's theorem makes it vanish for any
  /// closed molecule mesh, so it checks the outline and the cells).
  double area_mismatch = 0;
  bool valid = false;  ///< no overlap: all of the above within tolerance
};

OverlapReport check_overlap(const MoleculeMesh& mesh, const Layout& layout);

/// The Delaunay triangle of a vertex-tucking molecule on one of its sides,
/// and the 2D crease-pattern conditions on it (TVCG10 Sec. 4.2.2 and 7.2).
struct SideCondition {
  std::size_t vertex = 0;  ///< the side of VTM(vertex) at this slot
  std::size_t slot = 0;
  /// The angle of the Delaunay triangle opposite the side, and its base
  /// angles at the side's ends (counterclockwise order).
  double phi = 0;
  double beta0 = 0;
  double beta1 = 0;
  /// T9's angle between the side and a diagonal of the edge-tucking
  /// molecule; pi / 2 for a boundary molecule's side, whose Voronoi vertex
  /// must not leave the paper.
  double gamma = 0;
  bool t9 = true;   ///< phi <= gamma + pi / 2: the Voronoi vertex stays in the molecule
  bool t19 = true;  ///< min(beta0, beta1) >= -theta(v, j) / 2: the arc stays in the triangle
};

/// T9 and T19 on every side of positive length of every vertex-tucking
/// molecule, from the Delaunay triangulation of its corners (CGAL). For a
/// boundary molecule's side, T9 becomes phi <= pi / 2 and T19 does not apply.
std::vector<SideCondition> crease_conditions(const MoleculeMesh& mesh, const Layout& layout,
                                             const std::vector<double>& theta,
                                             const std::vector<double>& w);

/// Places the polygons from the angles and widths (TVCG10 Sec. 3 and 4.1):
/// the first face at its own coordinates, then every polygon across an edge
/// or digon molecule from a placed one, by the geometry of the symmetric
/// edge-tucking molecule (docs/tvcg2010-notes, Fig. 1). Then checks every
/// side of every vertex-tucking molecule, which tests T15 and T17
/// independently of how they were assembled. Throws InvariantError if the
/// closure error exceeds 1e-6 of the mean molecule length.
Layout reconstruct_layout(const SurfaceMesh& disk, const MoleculeMesh& mesh,
                          const std::vector<double>& theta, const std::vector<double>& w);

struct MappingOptions {
  /// Cuts to try in order (docs/notes.tex, M2): the fewest seams first, then
  /// every vertex on the cut.
  std::vector<CutStrategy> cuts = {CutStrategy::kMinimal, CutStrategy::kTree};
  std::size_t root_face = 0;
  MoleculeOptions molecules;
  RotationOptions rotation;  ///< its reference is replaced when remapping
  TranslationOptions translation;
};

/// One cut tried by map_surface, with or without remapping.
struct MappingAttempt {
  CutStrategy cut = CutStrategy::kMinimal;
  bool remapped = false;
  bool feasible = false;
  std::string message;  ///< why the attempt failed, from the step that failed
};

/// A surface mapped to the plane: the cut, the molecule mesh, the angles and
/// widths, and every attempt that led there.
struct Mapping {
  bool feasible = false;
  CutResult cut;
  MoleculeMesh molecules;
  std::vector<double> theta;  ///< from the last attempt that got past the rotation step
  std::vector<double> w;      ///< its widths, negative ones included if infeasible
  std::vector<MappingAttempt> attempts;
  /// The placed polygons, when the mapping is feasible; map_surface checks
  /// that they close and do not overlap.
  Layout layout;
};

/// The two-step mapping with its fallbacks (TVCG10 Sec. 6, Algorithm 1, as
/// adapted in docs/notes.tex, M2). For each cut in turn: build the molecule
/// mesh, run the rotation step and the translation step, and if the
/// translation fails, remap (remap_angles) and run both steps again. Stops at
/// the first attempt with nonnegative widths. If none succeeds, the result
/// holds the last attempt that reached the translation step, with phase 1's
/// widths when objective B is used, for diagnosis.
Mapping map_surface(const SurfaceMesh& mesh, const MappingOptions& options = {});

/// Numerical rank of the matrix of a set of equations (column-pivoted QR, with
/// a threshold of 1e-10 relative to the largest pivot). Used for the
/// degree-of-freedom check of TVCG10 Sec. 6.1 (docs/tvcg2010-notes, Open
/// question 5).
std::size_t numerical_rank(const LinearEquations& equations);

}  // namespace origamizer
