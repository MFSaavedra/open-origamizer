// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

#include "origamizer/layout.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <string>

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

namespace {

/// The elastic version of an infeasible rotation program (docs/notes.tex,
/// M2): the equality rows (T15) stay; every variable bound and inequality row
/// (T4, T6, T8) gets a slack s >= 0 that widens it on both sides, and the
/// total slack is minimized, a linear program. Fills in the least-violation
/// angles and the broken inequalities. Rows after the T15 rows are the T8
/// rows, star by star and slot by slot, as solve_rotation adds them.
void least_violation(const MoleculeMesh& mesh, const QuadraticProgram& program,
                     RotationResult& result) {
  const std::size_t n = program.variables();
  const std::size_t equalities = mesh.stars.size();
  QuadraticProgram elastic;
  for (std::size_t k = 0; k < n; ++k) elastic.add_variable(-kInfinity, kInfinity, 0.0);
  for (std::size_t r = 0; r < equalities; ++r) {
    elastic.add_row(program.row_lower[r], program.row_upper[r]);
  }
  for (const auto& e : program.entries) {
    if (e.row < equalities) elastic.entries.push_back(e);
  }
  // A two-sided constraint lo <= q <= hi becomes q + s >= lo and q - s <= hi.
  const auto relax = [&](const std::vector<QuadraticProgram::Entry>& terms, double lo, double hi) {
    const std::size_t slack = elastic.add_variable(0.0, kInfinity, 1.0);
    if (std::isfinite(lo)) {
      const std::size_t row = elastic.add_row(lo, kInfinity);
      for (const auto& t : terms) elastic.entries.push_back({row, t.column, t.value});
      elastic.entries.push_back({row, slack, 1.0});
    }
    if (std::isfinite(hi)) {
      const std::size_t row = elastic.add_row(-kInfinity, hi);
      for (const auto& t : terms) elastic.entries.push_back({row, t.column, t.value});
      elastic.entries.push_back({row, slack, -1.0});
    }
    return slack;
  };
  std::vector<std::size_t> bound_slack(n), row_slack(program.rows(), 0);
  for (std::size_t k = 0; k < n; ++k) {
    bound_slack[k] = relax({{0, k, 1.0}}, program.lower[k], program.upper[k]);
  }
  std::vector<std::vector<QuadraticProgram::Entry>> rows(program.rows());
  for (const auto& e : program.entries) rows[e.row].push_back(e);
  for (std::size_t r = equalities; r < program.rows(); ++r) {
    row_slack[r] = relax(rows[r], program.row_lower[r], program.row_upper[r]);
  }

  const detail::QpSolution solution = detail::solve(elastic);
  if (solution.status != detail::SolveStatus::kOptimal) {
    result.message += "; the elastic program failed too: " + solution.message;
    return;
  }
  result.theta.assign(solution.x.begin(), solution.x.begin() + static_cast<std::ptrdiff_t>(n));
  constexpr double kReported = 1e-9;
  for (std::size_t k = 0; k < n; ++k) {
    const double amount = solution.x[bound_slack[k]];
    if (amount <= kReported) continue;
    RotationViolation violation;
    violation.kind = mesh.molecules[k].kind == Molecule::Kind::kBoundary
                         ? RotationViolation::Kind::kT4
                         : RotationViolation::Kind::kT6;
    violation.molecule = k;
    violation.amount = amount;
    result.violations.push_back(violation);
  }
  std::size_t r = equalities;
  for (std::size_t v = 0; v < mesh.stars.size(); ++v) {
    for (std::size_t slot = 0; slot < mesh.stars[v].slots.size(); ++slot, ++r) {
      const double amount = solution.x[row_slack[r]];
      if (amount <= kReported) continue;
      RotationViolation violation;
      violation.kind = RotationViolation::Kind::kT8;
      violation.vertex = v;
      violation.slot = slot;
      violation.amount = amount;
      result.violations.push_back(violation);
    }
  }
  ORIGAMIZER_CHECK(r == program.rows(), "the T8 rows do not match the stars");
  std::sort(
      result.violations.begin(), result.violations.end(),
      [](const RotationViolation& a, const RotationViolation& b) { return a.amount > b.amount; });
  for (const RotationViolation& violation : result.violations) result.violation += violation.amount;
  result.message += "; least total violation " + format(result.violation) + " rad in " +
                    std::to_string(result.violations.size()) + " inequalities";
  for (std::size_t k = 0; k < std::min<std::size_t>(3, result.violations.size()); ++k) {
    const RotationViolation& violation = result.violations[k];
    result.message += k == 0 ? ", largest: " : ", ";
    if (violation.kind == RotationViolation::Kind::kT8) {
      result.message += "T8 at vertex " + std::to_string(violation.vertex) + " slot " +
                        std::to_string(violation.slot);
    } else {
      result.message += std::string(violation.kind == RotationViolation::Kind::kT4 ? "T4" : "T6") +
                        " at molecule " + std::to_string(violation.molecule);
    }
    result.message += " by " + format(violation.amount);
  }
}

}  // namespace

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
  if (solution.status == detail::SolveStatus::kInfeasible) {
    least_violation(mesh, program, result);
    return result;
  }
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

namespace {

Eigen::MatrixXd dense_matrix(const LinearEquations& equations) {
  Eigen::MatrixXd a = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(equations.rhs.size()),
                                            static_cast<Eigen::Index>(equations.columns));
  for (const SparseEntry& e : equations.entries) {
    ORIGAMIZER_CHECK(e.row < equations.rhs.size() && e.column < equations.columns,
                     "a matrix entry is out of range");
    a(static_cast<Eigen::Index>(e.row), static_cast<Eigen::Index>(e.column)) += e.value;
  }
  return a;
}

/// Adds the rows of `equations` as equality rows of `program`, on its
/// variables 0 to equations.columns - 1.
void add_equations(QuadraticProgram& program, const LinearEquations& equations) {
  const std::size_t first = program.rows();
  for (const double rhs : equations.rhs) program.add_row(rhs, rhs);
  for (const SparseEntry& e : equations.entries) {
    program.entries.push_back({first + e.row, e.column, e.value});
  }
}

std::vector<double> residual_vector(const LinearEquations& equations,
                                    const std::vector<double>& x) {
  std::vector<double> r(equations.rhs.size(), 0);
  for (const SparseEntry& e : equations.entries) r[e.row] += e.value * x[e.column];
  for (std::size_t k = 0; k < r.size(); ++k) r[k] -= equations.rhs[k];
  return r;
}

/// A maximal linearly independent subset of the rows of `equations` (the
/// first rank columns of a column-pivoted QR of the transpose), renumbered.
/// The active-set QP solver of HiGHS needs independent equality rows.
LinearEquations independent_rows(const LinearEquations& equations, const Eigen::MatrixXd& a) {
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(a.transpose());
  qr.setThreshold(1e-10);
  constexpr std::size_t kDropped = static_cast<std::size_t>(-1);
  std::vector<std::size_t> renumber(equations.rhs.size(), kDropped);
  LinearEquations independent;
  independent.columns = equations.columns;
  for (Eigen::Index k = 0; k < qr.rank(); ++k) {
    const auto row = static_cast<std::size_t>(qr.colsPermutation().indices()(k));
    renumber[row] = independent.rhs.size();
    independent.rhs.push_back(equations.rhs[row]);
  }
  for (const SparseEntry& e : equations.entries) {
    if (renumber[e.row] != kDropped)
      independent.entries.push_back({renumber[e.row], e.column, e.value});
  }
  return independent;
}

}  // namespace

double mean_molecule_length(const MoleculeMesh& mesh) {
  double sum = 0;
  std::size_t count = 0;
  for (const Molecule& m : mesh.molecules) {
    if (m.kind == Molecule::Kind::kBoundary) continue;
    sum += m.length;
    ++count;
  }
  ORIGAMIZER_CHECK(count > 0, "the molecule mesh has no edge or digon molecules");
  return sum / static_cast<double>(count);
}

std::vector<std::vector<double>> width_excess(const MoleculeMesh& mesh,
                                              const std::vector<double>& theta) {
  ORIGAMIZER_CHECK(theta.size() == mesh.molecules.size(), "expected one angle per molecule");
  std::vector<std::vector<double>> excess(mesh.stars.size());
  for (std::size_t v = 0; v < mesh.stars.size(); ++v) {
    for (const StarSlot& slot : mesh.stars[v].slots) {
      const Molecule& m = mesh.molecules[slot.molecule];
      const double angle = slot.sign * theta[slot.molecule];
      excess[v].push_back(m.kind == Molecule::Kind::kBoundary
                              ? 0.0
                              : 2 * m.length * std::max(0.0, -std::sin(angle / 2)));
    }
  }
  return excess;
}

LinearEquations translation_equations(const MoleculeMesh& mesh, const std::vector<double>& theta) {
  const auto exterior = exterior_angles(mesh, theta);
  const auto excess = width_excess(mesh, theta);
  LinearEquations equations;
  equations.columns = mesh.molecules.size();
  equations.rhs.assign(2 * mesh.stars.size(), 0);
  for (std::size_t v = 0; v < mesh.stars.size(); ++v) {
    double direction = 0;
    const auto& slots = mesh.stars[v].slots;
    for (std::size_t n = 0; n < slots.size(); ++n) {
      if (n > 0) direction += exterior[v][n];
      const double c = std::cos(direction), s = std::sin(direction);
      equations.entries.push_back({2 * v, slots[n].molecule, c});
      equations.entries.push_back({2 * v + 1, slots[n].molecule, s});
      equations.rhs[2 * v] -= excess[v][n] * c;
      equations.rhs[2 * v + 1] -= excess[v][n] * s;
    }
  }
  return equations;
}

TranslationResiduals translation_residuals(const MoleculeMesh& mesh,
                                           const std::vector<double>& theta,
                                           const std::vector<double>& w) {
  ORIGAMIZER_CHECK(w.size() == mesh.molecules.size(), "expected one width per molecule");
  TranslationResiduals residuals;
  for (const double r : residual_vector(translation_equations(mesh, theta), w)) {
    residuals.equations = std::max(residuals.equations, std::abs(r));
  }
  residuals.smallest_width = *std::min_element(w.begin(), w.end());
  return residuals;
}

TranslationResult solve_translation(const MoleculeMesh& mesh, const std::vector<double>& theta,
                                    const TranslationOptions& options) {
  const std::size_t n = mesh.molecules.size();
  // Solve in units of the mean molecule length, so that the solver's absolute
  // tolerances (1e-7) stay small relative to the widths on meshes of any size.
  const double scale = mean_molecule_length(mesh);
  LinearEquations all_rows = translation_equations(mesh, theta);
  for (double& rhs : all_rows.rhs) rhs /= scale;
  constexpr double kTolerance = 1e-6;  // relative to the scale
  TranslationResult result;

  // T17 alone: its minimum-norm least-squares solution (the paper's starting
  // point) and whether it solves the equations. On some meshes the rows are
  // dependent but consistent (docs/notes.tex, M2); the programs get an
  // independent subset, and the final check covers every row.
  const Eigen::MatrixXd a = dense_matrix(all_rows);
  const Eigen::VectorXd b = Eigen::Map<const Eigen::VectorXd>(
      all_rows.rhs.data(), static_cast<Eigen::Index>(all_rows.rhs.size()));
  const Eigen::VectorXd minimum_norm = a.completeOrthogonalDecomposition().solve(b);
  if (const double residual = (a * minimum_norm - b).lpNorm<Eigen::Infinity>();
      residual > kTolerance) {
    result.message =
        "T17 has no solution (least-squares residual " + format(residual * scale) + ")";
    return result;
  }
  const LinearEquations equations = independent_rows(all_rows, a);

  // Converts scaled widths back, checks them against the unscaled T17 and the
  // (scaled) lower bound, and stores them.
  const auto accept = [&](const std::vector<double>& x, double lower) {
    std::vector<double> w(n);
    for (std::size_t k = 0; k < n; ++k) w[k] = x[k] * scale;
    const TranslationResiduals r = translation_residuals(mesh, theta, w);
    ORIGAMIZER_CHECK(
        r.equations < kTolerance * scale && r.smallest_width > (lower - kTolerance) * scale,
        "the solver's widths violate the translation constraints (T17 by " + format(r.equations) +
            ", smallest width " + format(r.smallest_width) + ")");
    result.w = std::move(w);
    result.feasible = true;
  };

  if (options.objective == WidthObjective::kMinimumWidth) {
    ORIGAMIZER_CHECK(std::isfinite(options.minimum_width), "the minimum width must be finite");
    const double lower = options.minimum_width < 0 ? 0.01 : options.minimum_width / scale;
    QuadraticProgram program;
    for (std::size_t k = 0; k < n; ++k) program.add_variable(lower, kInfinity, 1.0);
    add_equations(program, equations);
    const detail::QpSolution solution = detail::solve(program);
    result.message = "objective A: " + solution.message;
    if (solution.status == detail::SolveStatus::kOptimal) accept(solution.x, lower);
    return result;
  }

  // Objective B, phase 1, a linear program: free widths w and slacks s >= 0
  // with w + s >= 0; minimize the sum of s, the total negative width. It is
  // zero exactly when the paper's E_(5,7) can vanish. (A sum of squares here
  // left HiGHS's QP solver with free variables of zero curvature, on which it
  // failed; docs/notes.tex, M2.)
  QuadraticProgram phase1;
  for (std::size_t k = 0; k < n; ++k) phase1.add_variable(-kInfinity, kInfinity, 0.0);
  for (std::size_t k = 0; k < n; ++k) phase1.add_variable(0.0, kInfinity, 1.0);
  add_equations(phase1, equations);
  for (std::size_t k = 0; k < n; ++k) {
    const std::size_t row = phase1.add_row(0.0, kInfinity);
    phase1.entries.push_back({row, k, 1.0});
    phase1.entries.push_back({row, n + k, 1.0});
  }
  const detail::QpSolution first = detail::solve(phase1);
  result.message = "objective B, phase 1: " + first.message;
  if (first.status != detail::SolveStatus::kOptimal) return result;
  double negative = 0;
  for (std::size_t k = n; k < 2 * n; ++k) negative += first.x[k];
  if (negative > kTolerance) {
    result.violation = negative * scale;
    result.w.resize(n);
    for (std::size_t k = 0; k < n; ++k) result.w[k] = first.x[k] * scale;
    result.message +=
        "; T5 and T7 cannot hold (total negative width " + format(result.violation) + ")";
    return result;
  }

  // Phase 2: the nonnegative widths closest to the reference.
  std::vector<double> reference(n);
  if (options.reference.empty()) {
    reference.assign(minimum_norm.data(), minimum_norm.data() + minimum_norm.size());
  } else {
    ORIGAMIZER_CHECK(options.reference.size() == n, "the reference needs one width per molecule");
    for (std::size_t k = 0; k < n; ++k) reference[k] = options.reference[k] / scale;
  }
  QuadraticProgram phase2;
  phase2.diagonal.assign(n, 2.0);  // (1/2) * 2 * w^2 - 2 ref w  =  |w - ref|^2 - |ref|^2
  for (std::size_t k = 0; k < n; ++k) {
    ORIGAMIZER_CHECK(std::isfinite(reference[k]), "the reference widths must be finite");
    phase2.cost.push_back(-2.0 * reference[k]);
    phase2.lower.push_back(0.0);
    phase2.upper.push_back(kInfinity);
  }
  add_equations(phase2, equations);
  const detail::QpSolution second = detail::solve(phase2);
  result.message = "objective B, phase 2: " + second.message;
  if (second.status == detail::SolveStatus::kOptimal) accept(second.x, 0.0);
  return result;
}

namespace {

/// Plane vectors for the reconstruction; converted to Point2 at the end.
using Vec = std::array<double, 2>;

Vec operator+(const Vec& a, const Vec& b) { return {a[0] + b[0], a[1] + b[1]}; }
Vec operator-(const Vec& a, const Vec& b) { return {a[0] - b[0], a[1] - b[1]}; }
Vec operator*(double s, const Vec& a) { return {s * a[0], s * a[1]}; }
double norm(const Vec& a) { return std::hypot(a[0], a[1]); }
Vec rotated(const Vec& a, double angle) {
  const double c = std::cos(angle), s = std::sin(angle);
  return {c * a[0] - s * a[1], s * a[0] + c * a[1]};
}
Vec direction(double angle) { return {std::cos(angle), std::sin(angle)}; }

/// A rigid motion of the plane: rotation by `angle`, then translation.
struct Placement {
  double angle = 0;
  Vec shift{0, 0};
  Vec operator()(const Vec& q) const { return rotated(q, angle) + shift; }
};

/// The rigid motion taking local points qa, qb to targets a, b (|qb - qa| =
/// |b - a| up to rounding).
Placement placement_from(const Vec& qa, const Vec& qb, const Vec& a, const Vec& b) {
  const Vec dq = qb - qa, d = b - a;
  Placement p;
  p.angle = std::atan2(d[1], d[0]) - std::atan2(dq[1], dq[0]);
  p.shift = a - rotated(qa, p.angle);
  return p;
}

}  // namespace

Layout reconstruct_layout(const SurfaceMesh& disk, const MoleculeMesh& mesh,
                          const std::vector<double>& theta, const std::vector<double>& w) {
  ORIGAMIZER_CHECK(theta.size() == mesh.molecules.size() && w.size() == mesh.molecules.size(),
                   "expected one angle and one width per molecule");
  ORIGAMIZER_CHECK(
      disk.number_of_faces() == mesh.faces && disk.number_of_vertices() == mesh.stars.size(),
      "the molecule mesh does not belong to this disk");
  const std::size_t polygons = mesh.polygons();
  constexpr std::size_t kNone = static_cast<std::size_t>(-1);

  // Local coordinates: each face in its own plane (counterclockwise seen from
  // the clean side), each digon as a segment along the x axis.
  Layout layout;
  layout.vertices.resize(polygons);
  std::vector<std::vector<Vec>> local(polygons);
  for (const auto f : disk.faces()) {
    const auto p = static_cast<std::size_t>(f);
    std::vector<Point3> points;
    for (const auto v : CGAL::vertices_around_face(disk.halfedge(f), disk)) {
      layout.vertices[p].push_back(static_cast<std::size_t>(v));
      points.push_back(disk.point(v));
    }
    Vector3 normal(0, 0, 0);  // Newell's method
    for (std::size_t k = 0; k < points.size(); ++k) {
      normal = normal + CGAL::cross_product(points[k] - CGAL::ORIGIN,
                                            points[(k + 1) % points.size()] - CGAL::ORIGIN);
    }
    const Vector3 x_axis =
        (points[1] - points[0]) / std::sqrt((points[1] - points[0]).squared_length());
    const Vector3 n = normal / std::sqrt(normal.squared_length());
    const Vector3 y_axis = CGAL::cross_product(n, x_axis);
    for (const Point3& q : points)
      local[p].push_back({(q - points[0]) * x_axis, (q - points[0]) * y_axis});
  }
  for (std::size_t k = 0; k < mesh.digon_molecules; ++k) {
    const Molecule& m = mesh.molecules[mesh.edge_molecules + k];
    layout.vertices[mesh.faces + k] = {m.i, m.j};
    local[mesh.faces + k] = {{0, 0}, {m.length, 0}};
  }
  const auto local_of = [&](std::size_t p, std::size_t v) -> const Vec& {
    for (std::size_t k = 0; k < layout.vertices[p].size(); ++k) {
      if (layout.vertices[p][k] == v) return local[p][k];
    }
    ORIGAMIZER_CHECK(
        false, "vertex " + std::to_string(v) + " is not a corner of polygon " + std::to_string(p));
    return local[p][0];  // unreachable
  };

  // Each edge or digon molecule joins the polygon before it (P0) and after it
  // (P1), counterclockwise around its vertex i.
  std::vector<std::array<std::size_t, 2>> sides(mesh.molecules.size(), {kNone, kNone});
  std::vector<std::vector<std::size_t>> adjacent(polygons);
  for (std::size_t v = 0; v < mesh.stars.size(); ++v) {
    const auto& slots = mesh.stars[v].slots;
    for (std::size_t n = 0; n < slots.size(); ++n) {
      const Molecule& m = mesh.molecules[slots[n].molecule];
      if (m.kind == Molecule::Kind::kBoundary || slots[n].sign != 1) continue;
      sides[slots[n].molecule] = {slots[n].polygon, slots[(n + 1) % slots.size()].polygon};
      adjacent[slots[n].polygon].push_back(slots[n].molecule);
      adjacent[slots[(n + 1) % slots.size()].polygon].push_back(slots[n].molecule);
    }
  }

  // Breadth-first placement across the molecules (docs/tvcg2010-notes,
  // Fig. 1): with u0 the direction of the leg i -> j in P0,
  //   A1 = A0 + w(i, j) R(theta / 2 + pi / 2) u0,   u1 = R(theta) u0.
  std::vector<Placement> placed(polygons);
  std::vector<bool> done(polygons, false);
  std::vector<std::size_t> queue = {0};
  done[0] = true;
  const auto excess = [&](std::size_t molecule, double angle) {
    return 2 * mesh.molecules[molecule].length * std::max(0.0, -std::sin(angle / 2));
  };
  for (std::size_t head = 0; head < queue.size(); ++head) {
    const std::size_t p = queue[head];
    for (const std::size_t molecule : adjacent[p]) {
      const Molecule& m = mesh.molecules[molecule];
      const auto [p0, p1] = sides[molecule];
      const std::size_t q = p == p0 ? p1 : p0;
      if (done[q]) continue;
      const double t = theta[molecule];
      const double wi = w[molecule] + excess(molecule, t);  // w(i, j)
      Vec a0, b0, a1, b1;
      if (p == p0) {
        a0 = placed[p](local_of(p, m.i));
        b0 = placed[p](local_of(p, m.j));
        const Vec u0 = (1 / norm(b0 - a0)) * (b0 - a0);
        a1 = a0 + wi * rotated(u0, t / 2 + std::numbers::pi / 2);
        b1 = a1 + m.length * rotated(u0, t);
        placed[q] = placement_from(local_of(q, m.i), local_of(q, m.j), a1, b1);
      } else {
        a1 = placed[p](local_of(p, m.i));
        b1 = placed[p](local_of(p, m.j));
        const Vec u0 = rotated((1 / norm(b1 - a1)) * (b1 - a1), -t);
        a0 = a1 - wi * rotated(u0, t / 2 + std::numbers::pi / 2);
        b0 = a0 + m.length * u0;
        placed[q] = placement_from(local_of(q, m.i), local_of(q, m.j), a0, b0);
      }
      done[q] = true;
      queue.push_back(q);
    }
  }
  ORIGAMIZER_CHECK(queue.size() == polygons, "the molecule mesh is not connected");
  layout.corners.resize(polygons);
  for (std::size_t p = 0; p < polygons; ++p) {
    for (const Vec& q : local[p]) {
      const Vec c = placed[p](q);
      layout.corners[p].emplace_back(c[0], c[1]);
    }
  }

  // Every side of every VTM against Eq. 2 and 3: side n has length
  // w(v, j_n) and direction phi_v + psi_n, with phi_v read off one edge or
  // digon molecule of the star.
  const auto exterior = exterior_angles(mesh, theta);
  const auto widening = width_excess(mesh, theta);
  for (std::size_t v = 0; v < mesh.stars.size(); ++v) {
    const auto& slots = mesh.stars[v].slots;
    std::vector<double> psi(slots.size(), 0.0);
    for (std::size_t n = 1; n < slots.size(); ++n) psi[n] = psi[n - 1] + exterior[v][n];
    std::size_t first = kNone;
    for (std::size_t n = 0; n < slots.size() && first == kNone; ++n) {
      if (mesh.molecules[slots[n].molecule].kind != Molecule::Kind::kBoundary) first = n;
    }
    ORIGAMIZER_CHECK(first != kNone,
                     "vertex " + std::to_string(v) + " has only a boundary molecule");
    const Molecule& m = mesh.molecules[slots[first].molecule];
    const std::size_t before = slots[first].polygon;
    const Vec leg = placed[before](local_of(before, m.i == v ? m.j : m.i)) -
                    placed[before](local_of(before, v));
    const Vec side =
        rotated(leg, slots[first].sign * theta[slots[first].molecule] / 2 + std::numbers::pi / 2);
    const double phi = std::atan2(side[1], side[0]) - psi[first];
    for (std::size_t n = 0; n < slots.size(); ++n) {
      const std::size_t a = slots[n].polygon, b = slots[(n + 1) % slots.size()].polygon;
      const Vec actual = placed[b](local_of(b, v)) - placed[a](local_of(a, v));
      const double length = w[slots[n].molecule] + widening[v][n];
      const Vec expected = length * direction(phi + psi[n]);
      layout.closure_error = std::max(layout.closure_error, norm(actual - expected));
    }
  }
  const double scale = mean_molecule_length(mesh);
  ORIGAMIZER_CHECK(layout.closure_error <= 1e-6 * scale,
                   "the molecule mesh does not close (error " + format(layout.closure_error) +
                       "): theta and w must satisfy T15 and T17");
  return layout;
}

std::size_t numerical_rank(const LinearEquations& equations) {
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(dense_matrix(equations));
  qr.setThreshold(1e-10);
  return static_cast<std::size_t>(qr.rank());
}

}  // namespace origamizer
