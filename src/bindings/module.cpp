// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

// Python bindings for the C++ core: the extension module origamizer._core.

#include <CGAL/version.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <array>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "origamizer/check.hpp"
#include "origamizer/cut.hpp"
#include "origamizer/layout.hpp"
#include "origamizer/mesh.hpp"
#include "origamizer/molecules.hpp"

namespace nb = nanobind;
namespace og = origamizer;
using namespace nb::literals;

namespace {

using Xyz = std::array<double, 3>;

std::vector<Xyz> to_xyz(const std::vector<og::Point3>& points) {
  std::vector<Xyz> out;
  out.reserve(points.size());
  for (const auto& p : points) out.push_back({p.x(), p.y(), p.z()});
  return out;
}

Xyz to_xyz(const og::Vector3& v) { return {v.x(), v.y(), v.z()}; }

std::vector<og::Point3> to_points(const std::vector<Xyz>& xyz) {
  std::vector<og::Point3> out;
  out.reserve(xyz.size());
  for (const auto& [x, y, z] : xyz) out.emplace_back(x, y, z);
  return out;
}

std::vector<std::vector<std::array<double, 2>>> to_xy(
    const std::vector<std::vector<og::Point2>>& polygons) {
  std::vector<std::vector<std::array<double, 2>>> out(polygons.size());
  for (std::size_t p = 0; p < polygons.size(); ++p) {
    for (const auto& q : polygons[p]) out[p].push_back({q.x(), q.y()});
  }
  return out;
}

std::string counts(std::size_t vertices, std::size_t faces) {
  return std::to_string(vertices) + " vertices, " + std::to_string(faces) + " faces";
}

}  // namespace

NB_MODULE(_core, m) {
  m.doc() = "C++ core of Open Origamizer.";

  nb::exception<og::InvariantError>(m, "InvariantError", PyExc_RuntimeError);

  m.def(
      "cgal_version", [] { return std::string(CGAL_VERSION_STR); },
      "Version of the CGAL headers the core was compiled against.");

  m.def(
      "_fail_check_for_testing",
      [](const std::string& message) { ORIGAMIZER_CHECK(false, message); }, "message"_a,
      "Fail an invariant check on purpose. Test hook for the InvariantError translation.");

  // Input ------------------------------------------------------------------

  nb::class_<og::PolygonSoup>(m, "PolygonSoup",
                              "Vertices and polygon faces as read from a file, before any checks.")
      .def(
          "__init__",
          [](og::PolygonSoup* self, const std::vector<Xyz>& points,
             std::vector<std::vector<std::size_t>> faces) {
            new (self) og::PolygonSoup{to_points(points), std::move(faces)};
          },
          "points"_a, "faces"_a)
      .def_prop_rw(
          "points", [](const og::PolygonSoup& s) { return to_xyz(s.points); },
          [](og::PolygonSoup& s, const std::vector<Xyz>& xyz) { s.points = to_points(xyz); },
          "Vertex coordinates as [x, y, z] lists (a copy).")
      .def_rw("faces", &og::PolygonSoup::faces,
              "Faces as lists of 0-based vertex indices (a copy).")
      .def("__repr__", [](const og::PolygonSoup& s) {
        return "PolygonSoup(" + counts(s.points.size(), s.faces.size()) + ")";
      });

  m.def(
      "parse_obj",
      [](const std::string& text, const std::string& source) {
        std::istringstream in(text);
        return og::parse_obj(in, source);
      },
      "text"_a, "source"_a = "<input>",
      "Parse OBJ text. Raises InvariantError, naming the line, for malformed input.");
  m.def("read_obj", &og::read_obj, "path"_a,
        "Read the vertices and polygon faces of an OBJ file. Raises InvariantError for malformed "
        "content and RuntimeError if the file cannot be opened.");

  // Checks -----------------------------------------------------------------

  nb::class_<og::SoupReport>(m, "SoupReport",
                             "Problems that keep a polygon soup from being one connected, "
                             "consistently oriented 2-manifold.")
      .def_ro("degenerate_faces", &og::SoupReport::degenerate_faces)
      .def_ro("duplicate_faces", &og::SoupReport::duplicate_faces)
      .def_ro("unused_vertices", &og::SoupReport::unused_vertices)
      .def_ro("coincident_vertices", &og::SoupReport::coincident_vertices)
      .def_ro("nonmanifold_edges", &og::SoupReport::nonmanifold_edges)
      .def_ro("nonmanifold_vertices", &og::SoupReport::nonmanifold_vertices)
      .def_ro("components", &og::SoupReport::components)
      .def_ro("orientable", &og::SoupReport::orientable)
      .def_ro("consistently_oriented", &og::SoupReport::consistently_oriented)
      .def("ok", &og::SoupReport::ok)
      .def("ok_except_orientation", &og::SoupReport::ok_except_orientation)
      .def("describe", &og::SoupReport::describe);
  m.def("check_soup", &og::check_soup, "soup"_a, "Run every combinatorial check on a soup.");

  nb::class_<og::FaceShape>(m, "FaceShape", "Shape of one face in its best-fit plane.")
      .def_prop_ro("normal", [](const og::FaceShape& s) { return to_xyz(s.normal); })
      .def_ro("area", &og::FaceShape::area)
      .def_ro("diameter", &og::FaceShape::diameter)
      .def_ro("planarity_error", &og::FaceShape::planarity_error)
      .def_ro("simple", &og::FaceShape::simple)
      .def_ro("convex", &og::FaceShape::convex);
  m.def("face_shape", &og::face_shape, "soup"_a, "face"_a);

  // Mesh -------------------------------------------------------------------

  nb::class_<og::Topology>(m, "Topology", "Counts that determine the topology of a surface.")
      .def_ro("vertices", &og::Topology::vertices)
      .def_ro("edges", &og::Topology::edges)
      .def_ro("faces", &og::Topology::faces)
      .def_ro("boundary_loops", &og::Topology::boundary_loops)
      .def_ro("components", &og::Topology::components)
      .def_ro("euler_characteristic", &og::Topology::euler_characteristic)
      .def_ro("genus", &og::Topology::genus)
      .def_prop_ro("closed", &og::Topology::closed)
      .def_prop_ro("is_disk", &og::Topology::is_disk)
      .def("__repr__", [](const og::Topology& t) {
        return "Topology(V=" + std::to_string(t.vertices) + ", E=" + std::to_string(t.edges) +
               ", F=" + std::to_string(t.faces) +
               ", chi=" + std::to_string(t.euler_characteristic) +
               ", boundary_loops=" + std::to_string(t.boundary_loops) +
               ", genus=" + std::to_string(t.genus) + ")";
      });

  nb::class_<og::SurfaceMesh>(m, "Mesh",
                              "A half-edge mesh with polygon faces (CGAL::Surface_mesh).")
      .def_prop_ro(
          "points", [](const og::SurfaceMesh& mesh) { return to_xyz(og::to_soup(mesh).points); },
          "Vertex coordinates as [x, y, z] lists.")
      .def_prop_ro(
          "faces", [](const og::SurfaceMesh& mesh) { return og::to_soup(mesh).faces; },
          "Faces as vertex index cycles, counterclockwise seen from the clean side.")
      .def("topology", &og::topology)
      .def("boundary_loops", &og::boundary_loops,
           "Boundary loops as vertex index cycles along the border halfedges.")
      .def(
          "gaussian_curvature",
          [](const og::SurfaceMesh& mesh) {
            std::vector<double> curvature;
            curvature.reserve(mesh.number_of_vertices());
            for (const auto v : mesh.vertices())
              curvature.push_back(og::gaussian_curvature(mesh, v));
            return curvature;
          },
          "Discrete Gaussian curvature of every vertex (angle defect inside, turning on the "
          "boundary).")
      .def("signed_volume", &og::signed_volume)
      .def("__repr__", [](const og::SurfaceMesh& mesh) {
        return "Mesh(" + counts(mesh.number_of_vertices(), mesh.number_of_faces()) + ")";
      });

  m.def("build_mesh", &og::build_mesh, "soup"_a, "allow_coincident_vertices"_a = false);
  m.def("to_soup", &og::to_soup, "mesh"_a);
  m.def("simplify_triangle_mesh", &og::simplify_triangle_mesh, "mesh"_a, "target_faces"_a,
        "Simplify a triangle mesh by edge collapses to at most target_faces faces, keeping its "
        "topology (CGAL, Lindstrom-Turk cost and placement).");

  nb::class_<og::PreparedMesh>(m, "PreparedMesh", "A checked mesh and what preparing it changed.")
      .def_ro("mesh", &og::PreparedMesh::mesh)
      .def_ro("face_origin", &og::PreparedMesh::face_origin,
              "For each mesh face, the face of the input it comes from.")
      .def_ro("split_faces", &og::PreparedMesh::split_faces,
              "Input faces that were split into convex pieces.")
      .def_ro("triangulated_faces", &og::PreparedMesh::triangulated_faces,
              "Nonplanar input quadrilaterals that were split into triangles.")
      .def_ro("reoriented_faces", &og::PreparedMesh::reoriented_faces,
              "Input faces reversed to make the orientation consistent.")
      .def_ro("flipped_to_clean_side", &og::PreparedMesh::flipped_to_clean_side,
              "Whether all faces were reversed so that a closed mesh faces outward.");

  m.def(
      "prepare_mesh",
      [](og::PolygonSoup soup, bool make_orientation_consistent, bool triangulate_nonplanar,
         double planarity_tolerance, double area_tolerance) {
        og::PrepareOptions options;
        options.make_orientation_consistent = make_orientation_consistent;
        options.triangulate_nonplanar = triangulate_nonplanar;
        options.tolerances.planarity = planarity_tolerance;
        options.tolerances.area = area_tolerance;
        return og::prepare_mesh(std::move(soup), options);
      },
      "soup"_a, nb::kw_only(), "make_orientation_consistent"_a = false,
      "triangulate_nonplanar"_a = false, "planarity_tolerance"_a = og::Tolerances{}.planarity,
      "area_tolerance"_a = og::Tolerances{}.area,
      "Check a soup, split nonconvex faces, build the mesh, and orient its clean side.");

  // Cutting ----------------------------------------------------------------

  nb::class_<og::CutResult>(m, "CutResult", "A surface cut into a topological disk.")
      .def_ro("disk", &og::CutResult::disk)
      .def_ro("original_vertex", &og::CutResult::original_vertex,
              "For each disk vertex, the mesh vertex it is a copy of.")
      .def_ro("cut_edges", &og::CutResult::cut_edges, "Mesh edges that were cut, as vertex pairs.");

  nb::enum_<og::CutStrategy>(m, "CutStrategy", "How much of the cut graph to keep.")
      .value("MINIMAL", og::CutStrategy::kMinimal, "the fewest cut edges")
      .value("TREE", og::CutStrategy::kTree,
             "the whole cut tree of TVCG10 Sec. 8.1: every vertex on the cut");

  m.def("automatic_cut", &og::automatic_cut, "mesh"_a, "root_face"_a = 0,
        "strategy"_a = og::CutStrategy::kMinimal,
        "Choose interior edges whose removal leaves a disk (TVCG10 Sec. 8.1).");
  m.def("cut_along", &og::cut_along, "mesh"_a, "edges"_a,
        "Cut the mesh along the given edges; the result must be a disk.");
  m.def(
      "cut_to_disk",
      [](const og::SurfaceMesh& mesh, std::size_t root_face, og::CutStrategy strategy,
         std::vector<og::IndexPair> cut_edges) {
        og::CutOptions options;
        options.root_face = root_face;
        options.strategy = strategy;
        options.cut_edges = std::move(cut_edges);
        return og::cut_to_disk(mesh, options);
      },
      "mesh"_a, nb::kw_only(), "root_face"_a = 0, "strategy"_a = og::CutStrategy::kMinimal,
      "cut_edges"_a = std::vector<og::IndexPair>{},
      "Cut along the given edges, or along the automatic cut when none are given.");

  // Molecule mesh (M2) -----------------------------------------------------

  nb::enum_<og::Molecule::Kind>(m, "MoleculeKind", "The three kinds of edge-tucking molecules.")
      .value("EDGE", og::Molecule::Kind::kEdge, "between the faces of an interior edge")
      .value("DIGON", og::Molecule::Kind::kDigon,
             "between a boundary face and the digon on its edge (TVCG10 Sec. 5.1)")
      .value("BOUNDARY", og::Molecule::Kind::kBoundary,
             "the outside of the paper at a boundary vertex");

  nb::class_<og::Molecule>(m, "Molecule", "An edge-tucking molecule, oriented from vertex i.")
      .def_ro("kind", &og::Molecule::kind)
      .def_ro("i", &og::Molecule::i)
      .def_ro("j", &og::Molecule::j, "The other vertex; equals i for a boundary molecule.")
      .def_ro("length", &og::Molecule::length)
      .def("__repr__", [](const og::Molecule& molecule) {
        const char* kind = molecule.kind == og::Molecule::Kind::kEdge    ? "edge"
                           : molecule.kind == og::Molecule::Kind::kDigon ? "digon"
                                                                         : "boundary";
        return std::string("Molecule(") + kind + ", " + std::to_string(molecule.i) + " -> " +
               std::to_string(molecule.j) + ")";
      });

  nb::class_<og::StarSlot>(m, "StarSlot", "One molecule in the star of a vertex.")
      .def_ro("molecule", &og::StarSlot::molecule)
      .def_ro("sign", &og::StarSlot::sign, "theta(v, j_n) = sign * theta of the molecule.")
      .def_ro("sector", &og::StarSlot::sector,
              "Corner angle of the polygon between the previous slot and this one.")
      .def_ro("polygon", &og::StarSlot::polygon,
              "That polygon: a face, or for a digon, faces plus its digon index.");

  nb::class_<og::VertexStar>(m, "VertexStar",
                             "The molecules around a vertex, counterclockwise from the clean side.")
      .def_ro("slots", &og::VertexStar::slots)
      .def_ro("boundary", &og::VertexStar::boundary)
      .def("defect", &og::VertexStar::defect);

  nb::class_<og::MoleculeMesh>(m, "MoleculeMesh", "The molecule mesh of a disk (TVCG10 Sec. 3).")
      .def_ro("molecules", &og::MoleculeMesh::molecules)
      .def_ro("stars", &og::MoleculeMesh::stars)
      .def_ro("edge_molecules", &og::MoleculeMesh::edge_molecules)
      .def_ro("digon_molecules", &og::MoleculeMesh::digon_molecules)
      .def_ro("faces", &og::MoleculeMesh::faces)
      .def_prop_ro("boundary_molecules", &og::MoleculeMesh::boundary_molecules)
      .def_prop_ro("polygons", &og::MoleculeMesh::polygons, "Faces, then digons.")
      .def("__repr__", [](const og::MoleculeMesh& mesh) {
        return "MoleculeMesh(" + std::to_string(mesh.molecules.size()) + " molecules, " +
               std::to_string(mesh.stars.size()) + " vertices)";
      });

  m.def(
      "build_molecule_mesh",
      [](const og::SurfaceMesh& disk, bool boundary_digons) {
        og::MoleculeOptions options;
        options.boundary_digons = boundary_digons;
        return og::build_molecule_mesh(disk, options);
      },
      "disk"_a, nb::kw_only(), "boundary_digons"_a = true,
      "Build the molecule mesh of a disk, with the boundary digons of TVCG10 Sec. 5.1 by "
      "default.");

  // Mapping (M2) -------------------------------------------------------------

  nb::enum_<og::RotationViolation::Kind>(m, "RotationViolationKind")
      .value("T4", og::RotationViolation::Kind::kT4, "theta(i, o) >= pi at a boundary molecule")
      .value("T6", og::RotationViolation::Kind::kT6, "|theta| <= pi - margin")
      .value("T8", og::RotationViolation::Kind::kT8, "0 <= Theta <= pi - margin at a corner");

  nb::class_<og::RotationViolation>(m, "RotationViolation",
                                    "An inequality that the least-violation angles break.")
      .def_ro("kind", &og::RotationViolation::kind)
      .def_ro("molecule", &og::RotationViolation::molecule, "For T4 and T6.")
      .def_ro("vertex", &og::RotationViolation::vertex, "For T8.")
      .def_ro("slot", &og::RotationViolation::slot, "For T8: the slot whose corner it is.")
      .def_ro("amount", &og::RotationViolation::amount, "In radians.");

  nb::class_<og::RotationResult>(m, "RotationResult", "The angles of the rotation step.")
      .def_ro("feasible", &og::RotationResult::feasible)
      .def_ro("theta", &og::RotationResult::theta,
              "One angle per molecule; when infeasible, the angles of least violation.")
      .def_ro("violation", &og::RotationResult::violation,
              "When infeasible: the least total violation of T4, T6, and T8 (radians).")
      .def_ro("violations", &og::RotationResult::violations, "The broken inequalities.")
      .def_ro("message", &og::RotationResult::message);

  m.def(
      "solve_rotation",
      [](const og::MoleculeMesh& mesh, double angle_margin, std::vector<double> reference) {
        og::RotationOptions options;
        options.angle_margin = angle_margin;
        options.reference = std::move(reference);
        return og::solve_rotation(mesh, options);
      },
      "mesh"_a, nb::kw_only(), "angle_margin"_a = og::RotationOptions{}.angle_margin,
      "reference"_a = std::vector<double>{},
      "The rotation step (TVCG10 Sec. 6.2): the angles closest to the reference (default zero) "
      "that satisfy T15, T4, T6, and T8.");
  m.def("remap_angles", &og::remap_angles, "mesh"_a,
        "Remapping (TVCG10 Sec. 6.5): reference angles from a Tutte drawing of the molecule "
        "graph.");
  m.def("exterior_angles", &og::exterior_angles, "mesh"_a, "theta"_a,
        "TVCG10 Eq. 3: the exterior angle of every vertex-tucking molecule at each slot.");
  m.def("mean_molecule_length", &og::mean_molecule_length, "mesh"_a);

  nb::enum_<og::WidthObjective>(m, "WidthObjective", "Objectives of the translation step.")
      .value("MINIMUM_WIDTH", og::WidthObjective::kMinimumWidth,
             "A: minimize the total width, each width at least a minimum")
      .value("PAPER_PENALTY", og::WidthObjective::kPaperPenalty,
             "B: TVCG10's penalty, as a feasibility phase and a closest-widths phase");

  nb::class_<og::TranslationResult>(m, "TranslationResult", "The widths of the translation step.")
      .def_ro("feasible", &og::TranslationResult::feasible)
      .def_ro("w", &og::TranslationResult::w,
              "One width per molecule; phase 1's widths when objective B is infeasible.")
      .def_ro("violation", &og::TranslationResult::violation,
              "Objective B: the smallest total negative width.")
      .def_ro("t19_failures", &og::TranslationResult::t19_failures,
              "With enforce_t19: the sides where T19 fails for the returned widths.")
      .def_ro("message", &og::TranslationResult::message);

  m.def(
      "solve_translation",
      [](const og::MoleculeMesh& mesh, const std::vector<double>& theta,
         og::WidthObjective objective, std::optional<double> minimum_width,
         std::vector<double> reference, bool enforce_t19) {
        og::TranslationOptions options;
        options.objective = objective;
        if (minimum_width) options.minimum_width = *minimum_width;
        options.reference = std::move(reference);
        options.enforce_t19 = enforce_t19;
        return og::solve_translation(mesh, theta, options);
      },
      "mesh"_a, "theta"_a, nb::kw_only(), "objective"_a = og::WidthObjective::kMinimumWidth,
      "minimum_width"_a = nb::none(), "reference"_a = std::vector<double>{},
      "enforce_t19"_a = false,
      "The translation step (TVCG10 Sec. 6.3) for given angles. minimum_width defaults to 1% of "
      "the mean molecule length; reference (objective B) to the minimum-norm solution of T17.");

  nb::class_<og::Layout>(m, "Layout", "The faces and digons of the molecule mesh in the plane.")
      .def_ro("vertices", &og::Layout::vertices,
              "For each polygon, the disk vertices at its corners.")
      .def_prop_ro(
          "corners", [](const og::Layout& layout) { return to_xy(layout.corners); },
          "For each polygon, its corners as [x, y] lists.")
      .def_ro("closure_error", &og::Layout::closure_error)
      .def_prop_ro(
          "outline", [](const og::Layout& layout) { return to_xy({layout.outline})[0]; },
          "The boundary of the paper, counterclockwise.")
      .def_prop_ro(
          "edge_molecules", [](const og::Layout& layout) { return to_xy(layout.edge_molecules); },
          "Each edge-tucking molecule as a quadrilateral, by molecule; empty for boundary "
          "molecules.")
      .def_prop_ro(
          "vertex_molecules",
          [](const og::Layout& layout) { return to_xy(layout.vertex_molecules); },
          "Each vertex-tucking molecule, its corners counterclockwise.");

  nb::class_<og::OverlapReport>(m, "OverlapReport",
                                "Whether a layout is a crease pattern without overlap.")
      .def_ro("smallest_area", &og::OverlapReport::smallest_area)
      .def_ro("degenerate_molecules", &og::OverlapReport::degenerate_molecules)
      .def_ro("reflex_turn", &og::OverlapReport::reflex_turn)
      .def_ro("outline_turning", &og::OverlapReport::outline_turning)
      .def_ro("area_mismatch", &og::OverlapReport::area_mismatch)
      .def_ro("valid", &og::OverlapReport::valid);
  nb::class_<og::SideCondition>(m, "SideCondition",
                                "T9 and T19 on one side of a vertex-tucking molecule.")
      .def_ro("vertex", &og::SideCondition::vertex)
      .def_ro("slot", &og::SideCondition::slot)
      .def_ro("phi", &og::SideCondition::phi)
      .def_ro("beta0", &og::SideCondition::beta0)
      .def_ro("beta1", &og::SideCondition::beta1)
      .def_ro("gamma", &og::SideCondition::gamma)
      .def_ro("t9", &og::SideCondition::t9)
      .def_ro("t19", &og::SideCondition::t19);
  m.def("crease_conditions", &og::crease_conditions, "mesh"_a, "layout"_a, "theta"_a, "w"_a,
        "T9 and T19 on every side of every vertex-tucking molecule (TVCG10 Sec. 4.2.2 and 7.2).");
  m.def("check_overlap", &og::check_overlap, "mesh"_a, "layout"_a,
        "Check that every molecule is convex with nonnegative area and the outline is convex "
        "(TVCG10 Sec. 4.2.1).");

  m.def("reconstruct_layout", &og::reconstruct_layout, "disk"_a, "mesh"_a, "theta"_a, "w"_a,
        "Place the polygons from the angles and widths and check that the molecule mesh closes.");

  nb::class_<og::MappingAttempt>(m, "MappingAttempt", "One cut tried by map_surface.")
      .def_ro("cut", &og::MappingAttempt::cut)
      .def_ro("remapped", &og::MappingAttempt::remapped)
      .def_ro("feasible", &og::MappingAttempt::feasible)
      .def_ro("message", &og::MappingAttempt::message)
      .def("__repr__", [](const og::MappingAttempt& a) {
        return std::string("MappingAttempt(") +
               (a.cut == og::CutStrategy::kMinimal ? "minimal" : "tree") + " cut" +
               (a.remapped ? ", remapped" : "") + (a.feasible ? ", feasible)" : ", infeasible)");
      });

  nb::class_<og::Mapping>(m, "Mapping", "A surface mapped to the plane, and how.")
      .def_ro("feasible", &og::Mapping::feasible)
      .def_ro("cut", &og::Mapping::cut)
      .def_ro("molecules", &og::Mapping::molecules)
      .def_ro("theta", &og::Mapping::theta)
      .def_ro("w", &og::Mapping::w)
      .def_ro("attempts", &og::Mapping::attempts)
      .def_ro("layout", &og::Mapping::layout, "The placed polygons, when feasible.");

  m.def(
      "map_surface",
      [](const og::SurfaceMesh& mesh, std::vector<og::CutStrategy> cuts, std::size_t root_face,
         bool boundary_digons, double angle_margin, og::WidthObjective objective,
         std::optional<double> minimum_width, bool enforce_t19) {
        og::MappingOptions options;
        options.cuts = std::move(cuts);
        options.root_face = root_face;
        options.molecules.boundary_digons = boundary_digons;
        options.rotation.angle_margin = angle_margin;
        options.translation.objective = objective;
        if (minimum_width) options.translation.minimum_width = *minimum_width;
        options.translation.enforce_t19 = enforce_t19;
        return og::map_surface(mesh, options);
      },
      "mesh"_a, nb::kw_only(), "cuts"_a = og::MappingOptions{}.cuts, "root_face"_a = 0,
      "boundary_digons"_a = true, "angle_margin"_a = og::RotationOptions{}.angle_margin,
      "objective"_a = og::WidthObjective::kMinimumWidth, "minimum_width"_a = nb::none(),
      "enforce_t19"_a = false,
      "The two-step mapping with its fallbacks: for each cut in turn, the rotation and "
      "translation steps, then remapping if the translation fails (TVCG10 Sec. 6).");
}
