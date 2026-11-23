// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

// Python bindings for the C++ core: the extension module origamizer._core.

#include <CGAL/version.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <array>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "origamizer/check.hpp"
#include "origamizer/cut.hpp"
#include "origamizer/mesh.hpp"

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

  nb::class_<og::PreparedMesh>(m, "PreparedMesh", "A checked mesh and what preparing it changed.")
      .def_ro("mesh", &og::PreparedMesh::mesh)
      .def_ro("face_origin", &og::PreparedMesh::face_origin,
              "For each mesh face, the face of the input it comes from.")
      .def_ro("split_faces", &og::PreparedMesh::split_faces,
              "Input faces that were split into convex pieces.")
      .def_ro("reoriented_faces", &og::PreparedMesh::reoriented_faces,
              "Input faces reversed to make the orientation consistent.")
      .def_ro("flipped_to_clean_side", &og::PreparedMesh::flipped_to_clean_side,
              "Whether all faces were reversed so that a closed mesh faces outward.");

  m.def(
      "prepare_mesh",
      [](og::PolygonSoup soup, bool make_orientation_consistent, double planarity_tolerance,
         double area_tolerance) {
        og::PrepareOptions options;
        options.make_orientation_consistent = make_orientation_consistent;
        options.tolerances.planarity = planarity_tolerance;
        options.tolerances.area = area_tolerance;
        return og::prepare_mesh(std::move(soup), options);
      },
      "soup"_a, nb::kw_only(), "make_orientation_consistent"_a = false,
      "planarity_tolerance"_a = og::Tolerances{}.planarity,
      "area_tolerance"_a = og::Tolerances{}.area,
      "Check a soup, split nonconvex faces, build the mesh, and orient its clean side.");

  // Cutting ----------------------------------------------------------------

  nb::class_<og::CutResult>(m, "CutResult", "A surface cut into a topological disk.")
      .def_ro("disk", &og::CutResult::disk)
      .def_ro("original_vertex", &og::CutResult::original_vertex,
              "For each disk vertex, the mesh vertex it is a copy of.")
      .def_ro("cut_edges", &og::CutResult::cut_edges, "Mesh edges that were cut, as vertex pairs.");

  m.def("automatic_cut", &og::automatic_cut, "mesh"_a, "root_face"_a = 0,
        "Choose interior edges whose removal leaves a disk (TVCG10 Sec. 8.1).");
  m.def("cut_along", &og::cut_along, "mesh"_a, "edges"_a,
        "Cut the mesh along the given edges; the result must be a disk.");
  m.def(
      "cut_to_disk",
      [](const og::SurfaceMesh& mesh, std::size_t root_face, std::vector<og::IndexPair> cut_edges) {
        og::CutOptions options;
        options.root_face = root_face;
        options.cut_edges = std::move(cut_edges);
        return og::cut_to_disk(mesh, options);
      },
      "mesh"_a, nb::kw_only(), "root_face"_a = 0, "cut_edges"_a = std::vector<og::IndexPair>{},
      "Cut along the given edges, or along the automatic cut when none are given.");
}
