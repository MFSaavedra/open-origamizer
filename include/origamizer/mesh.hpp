// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

/// \file
/// Mesh input and preparation (milestone M1).
///
/// The pipeline reads an OBJ file into a PolygonSoup, checks that it is one
/// connected, consistently oriented 2-manifold (check_soup), checks the shape
/// of every face, splits nonconvex faces into convex pieces, builds a
/// CGAL::Surface_mesh, and orients its clean side. prepare_mesh runs all of
/// these steps. Problems are reported, never repaired silently: the only
/// changes are the splitting of nonconvex faces, which the method requires,
/// and what the options explicitly ask for.

#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <cstddef>
#include <filesystem>
#include <iosfwd>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace origamizer {

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point2 = Kernel::Point_2;
using Vector2 = Kernel::Vector_2;
using Point3 = Kernel::Point_3;
using Vector3 = Kernel::Vector_3;
using SurfaceMesh = CGAL::Surface_mesh<Point3>;
using IndexPair = std::pair<std::size_t, std::size_t>;

// Input ----------------------------------------------------------------------

/// Vertices and polygon faces as read from a file, before any checks. Faces
/// hold 0-based vertex indices in their stored order, counterclockwise when
/// seen from the side the normal points to.
struct PolygonSoup {
  std::vector<Point3> points;
  std::vector<std::vector<std::size_t>> faces;
};

/// Parses Wavefront OBJ text. Reads `v` lines and `f` lines in the forms `i`,
/// `i/t`, `i/t/n`, and `i//n`, with negative (relative) indices allowed.
/// Texture coordinates, normals, object and group names, smoothing groups,
/// and materials are ignored. Any other element, such as `l` or `p`, is an
/// error, so that no geometry is dropped silently. Throws InvariantError
/// naming `source` and the line.
PolygonSoup parse_obj(std::istream& in, std::string_view source = "<input>");

/// Reads an OBJ file with parse_obj.
PolygonSoup read_obj(const std::filesystem::path& path);

// Combinatorial checks -------------------------------------------------------

/// Problems that keep a soup from being one connected, consistently oriented
/// 2-manifold, possibly with boundary. Indices refer to the soup.
struct SoupReport {
  /// Faces with fewer than 3 vertices, a repeated vertex, or an index out of range.
  std::vector<std::size_t> degenerate_faces;
  /// Pairs of faces on the same set of vertices.
  std::vector<IndexPair> duplicate_faces;
  /// Vertices that no face uses.
  std::vector<std::size_t> unused_vertices;
  /// Pairs of distinct vertices at exactly the same position.
  std::vector<IndexPair> coincident_vertices;
  /// Edges, as sorted vertex pairs, shared by more than two faces.
  std::vector<IndexPair> nonmanifold_edges;
  /// Vertices whose faces do not form a single fan.
  std::vector<std::size_t> nonmanifold_vertices;
  /// Connected components of the faces, joined through shared edges.
  std::size_t components = 0;
  /// Whether the faces can be oriented consistently (checked on manifold edges).
  bool orientable = true;
  /// Whether every edge shared by two faces is traversed in opposite directions.
  bool consistently_oriented = true;

  /// True when the only possible problem is an orientation that can be fixed.
  bool ok_except_orientation() const;
  /// True when the soup is one connected, consistently oriented 2-manifold.
  bool ok() const;
  /// Every problem, one per line, for error messages; empty when ok().
  std::string describe() const;
};

/// Runs every combinatorial check on the soup.
SoupReport check_soup(const PolygonSoup& soup);

/// Reverses faces so that every edge shared by two faces is traversed in
/// opposite directions, keeping the orientation of face 0. Requires a soup
/// whose only problem is a fixable orientation (ok_except_orientation() and
/// orientable). Returns the indices of the reversed faces.
std::vector<std::size_t> make_orientation_consistent(PolygonSoup& soup);

// Face shape -----------------------------------------------------------------

/// Tolerances of the geometric checks, relative to each face's diameter.
struct Tolerances {
  /// Largest allowed distance of a vertex from the face plane, over the diameter.
  double planarity = 1e-6;
  /// Smallest allowed face area, over the diameter squared.
  double area = 1e-12;
};

/// Shape of one face, measured in its best-fit (Newell) plane.
struct FaceShape {
  Vector3 normal = Vector3(0, 0, 0);  ///< unit normal, oriented by the vertex order
  double area = 0;
  double diameter = 0;         ///< largest distance between two of its vertices
  double planarity_error = 0;  ///< largest vertex distance from the plane, over the diameter
  bool simple = false;         ///< the projected boundary does not cross itself
  bool convex = false;         ///< simple with no reflex corner; straight corners allowed
};

/// Measures face `face` of the soup. Requires a non-degenerate face (see
/// SoupReport::degenerate_faces).
FaceShape face_shape(const PolygonSoup& soup, std::size_t face);

/// Result of split_nonconvex_faces.
struct SplitResult {
  PolygonSoup soup;                 ///< the same vertices; nonconvex faces replaced by their pieces
  std::vector<std::size_t> origin;  ///< face of `soup` -> face of the input
  std::vector<std::size_t> split;   ///< input faces that were split
};

/// Splits every nonconvex face into the fewest convex pieces, using diagonals
/// only (CGAL::optimal_convex_partition_2 in the face plane). The method needs
/// convex facets: concave corners violate the necessary conditions T13 and
/// T14 of TVCG10 Sec. 5.2, and SoCG17 Sec. 2 assumes convex facets. Throws
/// InvariantError listing every face that is too small, nonplanar, or not
/// simple under `tol`.
SplitResult split_nonconvex_faces(const PolygonSoup& soup, const Tolerances& tol = {});

// Mesh -----------------------------------------------------------------------

/// Builds a half-edge mesh from a soup that passes check_soup. Vertex i and
/// face i of the mesh are vertex i and face i of the soup. Throws
/// InvariantError listing every problem otherwise. Set
/// `allow_coincident_vertices` only for cut surfaces (cut.hpp), where the
/// copies of a split vertex share its position.
SurfaceMesh build_mesh(const PolygonSoup& soup, bool allow_coincident_vertices = false);

/// The vertex cycle of every face, as indices; the inverse of build_mesh up to
/// a rotation of each cycle.
PolygonSoup to_soup(const SurfaceMesh& mesh);

/// Counts that determine the topology of a surface.
struct Topology {
  std::size_t vertices = 0;
  std::size_t edges = 0;
  std::size_t faces = 0;
  std::size_t boundary_loops = 0;
  std::size_t components = 0;
  long euler_characteristic = 0;  ///< V - E + F
  long genus = 0;  ///< (2 c - chi - b) / 2: the genus of a connected orientable surface

  bool closed() const { return boundary_loops == 0; }
  bool is_disk() const {
    return components == 1 && boundary_loops == 1 && euler_characteristic == 1;
  }
};

/// Topology of an oriented mesh.
Topology topology(const SurfaceMesh& mesh);

/// Every boundary loop as a cycle of vertex indices. Each loop follows the
/// border halfedges, which run against the orientation of the faces; loops are
/// ordered by their first border halfedge.
std::vector<std::vector<std::size_t>> boundary_loops(const SurfaceMesh& mesh);

/// Sum of the face corner angles at vertex v.
double corner_angle_sum(const SurfaceMesh& mesh, SurfaceMesh::Vertex_index v);

/// Discrete Gaussian curvature at v: the angle defect 2 pi - sum(alpha) at an
/// interior vertex and the boundary turning pi - sum(alpha) at a boundary
/// vertex, so that the sum over all vertices is 2 pi chi (discrete
/// Gauss-Bonnet). TVCG10 Eq. 1 uses 2 pi - sum(alpha) at every vertex,
/// boundary vertices included (docs/tvcg2010-notes, Sec. 4).
double gaussian_curvature(const SurfaceMesh& mesh, SurfaceMesh::Vertex_index v);

/// Volume enclosed by a closed mesh; positive when the faces point outward.
double signed_volume(const SurfaceMesh& mesh);

/// Orients the faces so that their normals point to the clean side. A closed
/// mesh is oriented outward, which puts the tucks inside (CLAUDE.md, M1;
/// TVCG10 Sec. 8.2 and Fig. 20). A mesh with boundary keeps its orientation,
/// which is the user's choice of clean side. Returns true if it reversed the
/// faces. Throws InvariantError for a closed mesh that encloses no volume.
bool orient_clean_side_out(SurfaceMesh& mesh);

/// Simplifies a triangle mesh by edge collapses until it has at most
/// `target_faces` faces, using CGAL's Lindstrom-Turk cost and placement
/// (CGAL::Surface_mesh_simplification::edge_collapse), which preserves the
/// topology. Dense meshes must be simplified before folding (TVCG10 Sec. 10).
/// The result is renumbered. Throws InvariantError unless every face is a
/// triangle.
SurfaceMesh simplify_triangle_mesh(const SurfaceMesh& mesh, std::size_t target_faces);

// Pipeline -------------------------------------------------------------------

struct PrepareOptions {
  /// Reverse faces to make an orientable but inconsistently oriented soup
  /// consistent. Off by default, so that such input is reported as an error.
  bool make_orientation_consistent = false;
  Tolerances tolerances;
};

struct PreparedMesh {
  SurfaceMesh mesh;
  std::vector<std::size_t> face_origin;       ///< mesh face -> face of the input soup
  std::vector<std::size_t> split_faces;       ///< input faces split into convex pieces
  std::vector<std::size_t> reoriented_faces;  ///< input faces reversed for consistency
  bool flipped_to_clean_side = false;         ///< all faces reversed so a closed mesh faces out
};

/// Checks the soup, splits nonconvex faces, builds the mesh, and orients its
/// clean side, in that order. Throws InvariantError at the first stage that
/// fails, listing every problem that stage found.
PreparedMesh prepare_mesh(PolygonSoup soup, const PrepareOptions& options = {});

}  // namespace origamizer
