// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

#include "origamizer/mesh.hpp"

#include <CGAL/Partition_traits_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/partition_2.h>
#include <CGAL/property_map.h>

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iterator>
#include <list>
#include <map>
#include <numbers>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

#include "detail/disjoint_sets.hpp"
#include "origamizer/check.hpp"

namespace origamizer {
namespace {

using detail::DisjointSets;

// Formatting -----------------------------------------------------------------

std::string format_item(std::size_t i) { return std::to_string(i); }

std::string format_item(const IndexPair& p) {
  return "(" + std::to_string(p.first) + ", " + std::to_string(p.second) + ")";
}

/// At most `limit` items, then the number left out.
template <class T>
std::string format_list(const std::vector<T>& items, std::size_t limit = 12) {
  std::string out;
  for (std::size_t k = 0; k < items.size() && k < limit; ++k) {
    if (k > 0) out += ", ";
    out += format_item(items[k]);
  }
  if (items.size() > limit) out += ", ... (" + std::to_string(items.size() - limit) + " more)";
  return out;
}

std::string format_number(double x) {
  std::ostringstream out;
  out.precision(3);
  out << x;
  return out.str();
}

// OBJ parsing ----------------------------------------------------------------

bool parse_double(std::string_view text, double& value) {
  if (!text.empty() && text.front() == '+') text.remove_prefix(1);
  const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
  return error == std::errc() && end == text.data() + text.size();
}

// Edges of a soup -------------------------------------------------------------

/// One face's use of an edge; `forward` means it runs from the smaller vertex
/// index to the larger one.
struct Incidence {
  std::size_t face;
  bool forward;
};

using EdgeMap = std::map<IndexPair, std::vector<Incidence>>;

EdgeMap edge_incidences(const PolygonSoup& soup, const std::vector<char>& usable) {
  EdgeMap edges;
  for (std::size_t f = 0; f < soup.faces.size(); ++f) {
    if (!usable[f]) continue;
    const auto& face = soup.faces[f];
    for (std::size_t k = 0; k < face.size(); ++k) {
      const std::size_t a = face[k];
      const std::size_t b = face[(k + 1) % face.size()];
      edges[{std::min(a, b), std::max(a, b)}].push_back({f, a < b});
    }
  }
  return edges;
}

/// Neighbors of each face through edges shared by exactly two faces, with
/// whether the two faces traverse that edge in the same direction.
std::vector<std::vector<std::pair<std::size_t, bool>>> manifold_adjacency(const EdgeMap& edges,
                                                                          std::size_t faces) {
  std::vector<std::vector<std::pair<std::size_t, bool>>> adjacent(faces);
  for (const auto& [edge, uses] : edges) {
    if (uses.size() != 2) continue;
    const bool same = uses[0].forward == uses[1].forward;
    adjacent[uses[0].face].push_back({uses[1].face, same});
    adjacent[uses[1].face].push_back({uses[0].face, same});
  }
  return adjacent;
}

// Face geometry ----------------------------------------------------------------

Point3 face_centroid(const PolygonSoup& soup, const std::vector<std::size_t>& face) {
  double x = 0, y = 0, z = 0;
  for (const std::size_t i : face) {
    x += soup.points[i].x();
    y += soup.points[i].y();
    z += soup.points[i].z();
  }
  const auto n = static_cast<double>(face.size());
  return {x / n, y / n, z / n};
}

/// Newell vector (twice the vector area) of a face, computed relative to `c`.
Vector3 newell_vector(const PolygonSoup& soup, const std::vector<std::size_t>& face,
                      const Point3& c) {
  double x = 0, y = 0, z = 0;
  for (std::size_t k = 0; k < face.size(); ++k) {
    const Vector3 p = soup.points[face[k]] - c;
    const Vector3 q = soup.points[face[(k + 1) % face.size()]] - c;
    x += (p.y() - q.y()) * (p.z() + q.z());
    y += (p.z() - q.z()) * (p.x() + q.x());
    z += (p.x() - q.x()) * (p.y() + q.y());
  }
  return {x, y, z};
}

/// Projects a face onto the plane through `c` with unit normal `n`. The axes
/// (u, v, n) are right-handed, so a face that is counterclockwise around n
/// stays counterclockwise in the plane.
std::vector<Point2> project_face(const PolygonSoup& soup, const std::vector<std::size_t>& face,
                                 const Point3& c, const Vector3& n) {
  const double ax = std::abs(n.x()), ay = std::abs(n.y()), az = std::abs(n.z());
  const Vector3 axis = (ax <= ay && ax <= az) ? Vector3(1, 0, 0)
                       : (ay <= az)           ? Vector3(0, 1, 0)
                                              : Vector3(0, 0, 1);
  Vector3 u = axis - (axis * n) * n;
  u = u / std::sqrt(u.squared_length());
  const Vector3 v = CGAL::cross_product(n, u);
  std::vector<Point2> points;
  points.reserve(face.size());
  for (const std::size_t i : face) {
    const Vector3 d = soup.points[i] - c;
    points.emplace_back(d * u, d * v);
  }
  return points;
}

bool same_cycle(const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
  if (a.size() != b.size() || a.empty()) return false;
  const auto start = std::find(b.begin(), b.end(), a[0]);
  if (start == b.end()) return false;
  const auto offset = static_cast<std::size_t>(start - b.begin());
  for (std::size_t k = 0; k < a.size(); ++k) {
    if (a[k] != b[(offset + k) % b.size()]) return false;
  }
  return true;
}

double bounding_box_extent(const SurfaceMesh& mesh) {
  std::array<double, 3> lo{}, hi{};
  bool first = true;
  for (const auto v : mesh.vertices()) {
    const Point3& p = mesh.point(v);
    const std::array<double, 3> xyz{p.x(), p.y(), p.z()};
    for (int k = 0; k < 3; ++k) {
      lo[k] = first ? xyz[k] : std::min(lo[k], xyz[k]);
      hi[k] = first ? xyz[k] : std::max(hi[k], xyz[k]);
    }
    first = false;
  }
  return std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]});
}

}  // namespace

// Input ----------------------------------------------------------------------

PolygonSoup parse_obj(std::istream& in, std::string_view source) {
  PolygonSoup soup;
  std::vector<std::size_t> face_lines;
  std::string line;
  std::size_t line_number = 0;
  const auto fail = [&](const std::string& what) {
    throw InvariantError(std::string(source) + ":" + std::to_string(line_number) + ": " + what);
  };

  while (std::getline(in, line)) {
    ++line_number;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (const auto hash = line.find('#'); hash != std::string::npos) line.erase(hash);
    std::istringstream fields(line);
    std::string keyword;
    if (!(fields >> keyword)) continue;

    if (keyword == "v") {
      std::array<double, 3> xyz{};
      for (double& c : xyz) {
        std::string token;
        if (!(fields >> token) || !parse_double(token, c)) {
          fail("expected three coordinates after 'v'");
        }
        if (!std::isfinite(c)) fail("vertex coordinate '" + token + "' is not finite");
      }
      std::string token;
      double ignored = 0;
      while (fields >> token) {  // optional w or color components
        if (!parse_double(token, ignored)) fail("unexpected '" + token + "' after the coordinates");
      }
      soup.points.emplace_back(xyz[0], xyz[1], xyz[2]);
    } else if (keyword == "f") {
      std::vector<std::size_t> face;
      std::string token;
      while (fields >> token) {
        const std::string_view text = std::string_view(token).substr(0, token.find('/'));
        long long value = 0;
        const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
        if (error != std::errc() || end != text.data() + text.size() || value == 0) {
          fail("bad vertex reference '" + token + "'");
        }
        const auto count = static_cast<long long>(soup.points.size());
        const long long index = value > 0 ? value - 1 : count + value;
        if (index < 0) fail("vertex reference '" + token + "' points before the first vertex");
        face.push_back(static_cast<std::size_t>(index));
      }
      soup.faces.push_back(std::move(face));
      face_lines.push_back(line_number);
    } else if (keyword == "vt" || keyword == "vn" || keyword == "vp" || keyword == "o" ||
               keyword == "g" || keyword == "s" || keyword == "usemtl" || keyword == "mtllib") {
      continue;
    } else {
      fail("unsupported OBJ element '" + keyword + "'");
    }
  }

  for (std::size_t f = 0; f < soup.faces.size(); ++f) {
    for (const std::size_t i : soup.faces[f]) {
      if (i >= soup.points.size()) {
        line_number = face_lines[f];
        fail("face uses vertex " + std::to_string(i + 1) + ", but the file defines only " +
             std::to_string(soup.points.size()));
      }
    }
  }
  return soup;
}

PolygonSoup read_obj(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path.string());
  return parse_obj(in, path.string());
}

// Combinatorial checks -------------------------------------------------------

bool SoupReport::ok_except_orientation() const {
  return degenerate_faces.empty() && duplicate_faces.empty() && unused_vertices.empty() &&
         coincident_vertices.empty() && nonmanifold_edges.empty() && nonmanifold_vertices.empty() &&
         components == 1;
}

bool SoupReport::ok() const { return ok_except_orientation() && consistently_oriented; }

std::string SoupReport::describe() const {
  std::ostringstream out;
  const auto line = [&out](std::size_t count, const char* what, const std::string& items) {
    if (count > 0) out << count << ' ' << what << ": " << items << '\n';
  };
  line(degenerate_faces.size(),
       "degenerate faces (fewer than 3 vertices, a repeated vertex, or a missing vertex)",
       format_list(degenerate_faces));
  line(duplicate_faces.size(), "pairs of faces on the same vertices", format_list(duplicate_faces));
  line(unused_vertices.size(), "vertices used by no face", format_list(unused_vertices));
  line(coincident_vertices.size(), "pairs of distinct vertices at the same position",
       format_list(coincident_vertices));
  line(nonmanifold_edges.size(), "edges shared by more than two faces",
       format_list(nonmanifold_edges));
  line(nonmanifold_vertices.size(), "vertices whose faces do not form a single fan",
       format_list(nonmanifold_vertices));
  if (components != 1) out << components << " connected components (expected 1)\n";
  if (!orientable) {
    out << "the surface is not orientable\n";
  } else if (!consistently_oriented) {
    out << "the faces are not consistently oriented (the surface is orientable)\n";
  }
  return out.str();
}

SoupReport check_soup(const PolygonSoup& soup) {
  SoupReport report;
  const std::size_t n = soup.points.size();
  const std::size_t m = soup.faces.size();
  for (const Point3& p : soup.points) {
    ORIGAMIZER_CHECK(std::isfinite(p.x()) && std::isfinite(p.y()) && std::isfinite(p.z()),
                     "a vertex has a coordinate that is not finite");
  }

  // Degenerate faces and unused vertices.
  std::vector<char> usable(m, 1);
  std::vector<char> used(n, 0);
  for (std::size_t f = 0; f < m; ++f) {
    const auto& face = soup.faces[f];
    bool bad = face.size() < 3 ||
               std::any_of(face.begin(), face.end(), [n](std::size_t i) { return i >= n; });
    if (!bad) {
      std::vector<std::size_t> sorted(face);
      std::sort(sorted.begin(), sorted.end());
      bad = std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end();
    }
    if (bad) {
      report.degenerate_faces.push_back(f);
      usable[f] = 0;
      continue;
    }
    for (const std::size_t i : face) used[i] = 1;
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (!used[i]) report.unused_vertices.push_back(i);
  }

  // Coincident vertices: equal coordinates, compared exactly.
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), std::size_t{0});
  const auto key = [&soup](std::size_t i) {
    const Point3& p = soup.points[i];
    return std::make_tuple(p.x(), p.y(), p.z(), i);
  };
  std::sort(order.begin(), order.end(),
            [&key](std::size_t a, std::size_t b) { return key(a) < key(b); });
  for (std::size_t k = 1, start = 0; k < n; ++k) {
    if (soup.points[order[k]] == soup.points[order[start]]) {
      report.coincident_vertices.push_back({order[start], order[k]});
    } else {
      start = k;
    }
  }
  std::sort(report.coincident_vertices.begin(), report.coincident_vertices.end());

  // Duplicate faces.
  std::map<std::vector<std::size_t>, std::size_t> first_with_vertices;
  for (std::size_t f = 0; f < m; ++f) {
    if (!usable[f]) continue;
    std::vector<std::size_t> sorted(soup.faces[f]);
    std::sort(sorted.begin(), sorted.end());
    const auto [it, inserted] = first_with_vertices.emplace(std::move(sorted), f);
    if (!inserted) report.duplicate_faces.push_back({it->second, f});
  }

  // Edges: non-manifold edges, components, orientation.
  const EdgeMap edges = edge_incidences(soup, usable);
  DisjointSets faces_joined(m);
  for (const auto& [edge, uses] : edges) {
    if (uses.size() > 2) report.nonmanifold_edges.push_back(edge);
    for (std::size_t k = 1; k < uses.size(); ++k) faces_joined.unite(uses[0].face, uses[k].face);
    if (uses.size() == 2 && uses[0].forward == uses[1].forward) {
      report.consistently_oriented = false;
    }
  }
  std::set<std::size_t> roots;
  for (std::size_t f = 0; f < m; ++f) {
    if (usable[f]) roots.insert(faces_joined.find(f));
  }
  report.components = roots.size();

  // Orientability: two-color the faces by whether each must be reversed.
  const auto adjacent = manifold_adjacency(edges, m);
  std::vector<int> reverse(m, -1);
  for (std::size_t start = 0; start < m; ++start) {
    if (!usable[start] || reverse[start] != -1) continue;
    reverse[start] = 0;
    std::queue<std::size_t> pending;
    pending.push(start);
    while (!pending.empty()) {
      const std::size_t f = pending.front();
      pending.pop();
      for (const auto& [g, same] : adjacent[f]) {
        const int wanted = same ? 1 - reverse[f] : reverse[f];
        if (reverse[g] == -1) {
          reverse[g] = wanted;
          pending.push(g);
        } else if (reverse[g] != wanted) {
          report.orientable = false;
        }
      }
    }
  }

  // Non-manifold vertices: the link of a vertex (one segment per face corner)
  // must be a single path or cycle.
  std::vector<std::vector<IndexPair>> links(n);
  for (std::size_t f = 0; f < m; ++f) {
    if (!usable[f]) continue;
    const auto& face = soup.faces[f];
    const std::size_t s = face.size();
    for (std::size_t k = 0; k < s; ++k) {
      links[face[k]].push_back({face[(k + s - 1) % s], face[(k + 1) % s]});
    }
  }
  for (std::size_t v = 0; v < n; ++v) {
    if (links[v].empty()) continue;
    std::map<std::size_t, std::size_t> local;
    for (const auto& [a, b] : links[v]) {
      local.emplace(a, local.size());
      local.emplace(b, local.size());
    }
    std::vector<int> degree(local.size(), 0);
    DisjointSets joined(local.size());
    for (const auto& [a, b] : links[v]) {
      const std::size_t ia = local.at(a), ib = local.at(b);
      ++degree[ia];
      ++degree[ib];
      joined.unite(ia, ib);
    }
    std::set<std::size_t> link_roots;
    for (std::size_t k = 0; k < local.size(); ++k) link_roots.insert(joined.find(k));
    const bool fan = std::all_of(degree.begin(), degree.end(), [](int d) { return d <= 2; });
    if (!fan || link_roots.size() != 1) report.nonmanifold_vertices.push_back(v);
  }
  return report;
}

std::vector<std::size_t> make_orientation_consistent(PolygonSoup& soup) {
  const SoupReport report = check_soup(soup);
  ORIGAMIZER_CHECK(
      report.ok_except_orientation() && report.orientable,
      "only an orientable, connected 2-manifold can be reoriented:\n" + report.describe());
  const std::size_t m = soup.faces.size();
  const auto adjacent = manifold_adjacency(edge_incidences(soup, std::vector<char>(m, 1)), m);
  std::vector<int> reverse(m, -1);
  reverse[0] = 0;
  std::queue<std::size_t> pending;
  pending.push(0);
  while (!pending.empty()) {
    const std::size_t f = pending.front();
    pending.pop();
    for (const auto& [g, same] : adjacent[f]) {
      if (reverse[g] == -1) {
        reverse[g] = same ? 1 - reverse[f] : reverse[f];
        pending.push(g);
      }
    }
  }
  std::vector<std::size_t> reversed;
  for (std::size_t f = 0; f < m; ++f) {
    ORIGAMIZER_CHECK(reverse[f] != -1, "a face was not reached while reorienting");
    if (reverse[f] == 1) {
      std::reverse(soup.faces[f].begin(), soup.faces[f].end());
      reversed.push_back(f);
    }
  }
  ORIGAMIZER_CHECK(check_soup(soup).ok(), "reorienting left the faces inconsistent");
  return reversed;
}

// Face shape -----------------------------------------------------------------

FaceShape face_shape(const PolygonSoup& soup, std::size_t face_index) {
  ORIGAMIZER_CHECK(face_index < soup.faces.size(), "face index out of range");
  const auto& face = soup.faces[face_index];
  ORIGAMIZER_CHECK(face.size() >= 3, "a face needs at least 3 vertices");
  for (const std::size_t i : face) {
    ORIGAMIZER_CHECK(i < soup.points.size(), "a face uses a missing vertex");
  }

  FaceShape shape;
  for (std::size_t a = 0; a < face.size(); ++a) {
    for (std::size_t b = a + 1; b < face.size(); ++b) {
      shape.diameter =
          std::max(shape.diameter,
                   std::sqrt(CGAL::squared_distance(soup.points[face[a]], soup.points[face[b]])));
    }
  }
  const Point3 c = face_centroid(soup, face);
  const Vector3 newell = newell_vector(soup, face, c);
  const double length = std::sqrt(newell.squared_length());
  shape.area = length / 2;
  if (length == 0 || shape.diameter == 0) return shape;

  shape.normal = newell / length;
  for (const std::size_t i : face) {
    shape.planarity_error =
        std::max(shape.planarity_error, std::abs((soup.points[i] - c) * shape.normal));
  }
  shape.planarity_error /= shape.diameter;

  const std::vector<Point2> points = project_face(soup, face, c, shape.normal);
  shape.simple = CGAL::is_simple_2(points.begin(), points.end(), Kernel());
  if (shape.simple) {
    const std::size_t s = points.size();
    bool reflex = false;
    for (std::size_t k = 0; k < s; ++k) {
      if (CGAL::orientation(points[(k + s - 1) % s], points[k], points[(k + 1) % s]) ==
          CGAL::RIGHT_TURN) {
        reflex = true;
      }
    }
    shape.convex = !reflex;
  }
  return shape;
}

SplitResult split_nonconvex_faces(const PolygonSoup& soup, const Tolerances& tol) {
  using PartitionTraits =
      CGAL::Partition_traits_2<Kernel, CGAL::Pointer_property_map<Point2>::const_type>;

  SplitResult result;
  result.soup.points = soup.points;
  std::vector<std::string> problems;
  for (std::size_t f = 0; f < soup.faces.size(); ++f) {
    const auto& face = soup.faces[f];
    const FaceShape shape = face_shape(soup, f);
    const std::string name = "face " + std::to_string(f);
    if (!(shape.area > tol.area * shape.diameter * shape.diameter)) {
      problems.push_back(name + " has (nearly) zero area");
      continue;
    }
    if (shape.planarity_error > tol.planarity) {
      problems.push_back(
          name + " is not planar: a vertex lies " + format_number(shape.planarity_error) +
          " face diameters from its plane (tolerance " + format_number(tol.planarity) + ")");
      continue;
    }
    if (!shape.simple) {
      problems.push_back(name + " crosses itself");
      continue;
    }
    if (shape.convex) {
      result.soup.faces.push_back(face);
      result.origin.push_back(f);
      continue;
    }

    const std::vector<Point2> points =
        project_face(soup, face, face_centroid(soup, face), shape.normal);
    const PartitionTraits traits(CGAL::make_property_map(points));
    std::vector<std::size_t> corners(face.size());
    std::iota(corners.begin(), corners.end(), std::size_t{0});
    std::list<PartitionTraits::Polygon_2> pieces;
    CGAL::optimal_convex_partition_2(corners.begin(), corners.end(), std::back_inserter(pieces),
                                     traits);
    ORIGAMIZER_CHECK(pieces.size() >= 2, "a nonconvex face was partitioned into one piece");

    double piece_area = 0;
    for (const auto& piece : pieces) {
      std::vector<std::size_t> vertices;
      for (auto it = piece.vertices_begin(); it != piece.vertices_end(); ++it) {
        vertices.push_back(face[*it]);
      }
      ORIGAMIZER_CHECK(vertices.size() >= 3, "a convex piece has fewer than 3 vertices");
      piece_area +=
          std::sqrt(newell_vector(soup, vertices, face_centroid(soup, vertices)).squared_length()) /
          2;
      result.soup.faces.push_back(std::move(vertices));
      result.origin.push_back(f);
    }
    ORIGAMIZER_CHECK(std::abs(piece_area - shape.area) <= 1e-9 * shape.area,
                     "the convex pieces of " + name + " do not cover it exactly");
    result.split.push_back(f);
  }
  if (!problems.empty()) {
    std::string message = "faces that cannot be facets:";
    for (const auto& problem : problems) message += "\n" + problem;
    throw InvariantError(message);
  }
  return result;
}

// Mesh -----------------------------------------------------------------------

SurfaceMesh build_mesh(const PolygonSoup& soup, bool allow_coincident_vertices) {
  namespace PMP = CGAL::Polygon_mesh_processing;
  SoupReport report = check_soup(soup);
  if (allow_coincident_vertices) report.coincident_vertices.clear();
  if (!report.ok()) {
    throw InvariantError(
        "the polygon soup is not one connected, consistently oriented 2-manifold:\n" +
        report.describe());
  }
  ORIGAMIZER_CHECK(PMP::is_polygon_soup_a_polygon_mesh(soup.faces),
                   "CGAL rejects a polygon soup that passed check_soup");
  SurfaceMesh mesh;
  PMP::polygon_soup_to_polygon_mesh(soup.points, soup.faces, mesh);
  ORIGAMIZER_CHECK(mesh.number_of_vertices() == soup.points.size() &&
                       mesh.number_of_faces() == soup.faces.size(),
                   "building the mesh lost vertices or faces");
  for (std::size_t i = 0; i < soup.points.size(); ++i) {
    ORIGAMIZER_CHECK(mesh.point(SurfaceMesh::Vertex_index(
                         static_cast<SurfaceMesh::size_type>(i))) == soup.points[i],
                     "building the mesh renumbered the vertices");
  }
  const PolygonSoup rebuilt = to_soup(mesh);
  for (std::size_t f = 0; f < soup.faces.size(); ++f) {
    ORIGAMIZER_CHECK(same_cycle(soup.faces[f], rebuilt.faces[f]),
                     "building the mesh renumbered or reoriented a face");
  }
  ORIGAMIZER_CHECK(mesh.is_valid(false), "CGAL reports an invalid half-edge structure");
  return mesh;
}

PolygonSoup to_soup(const SurfaceMesh& mesh) {
  ORIGAMIZER_CHECK(!mesh.has_garbage(), "the mesh has removed elements; collect garbage first");
  PolygonSoup soup;
  soup.points.reserve(mesh.number_of_vertices());
  for (const auto v : mesh.vertices()) soup.points.push_back(mesh.point(v));
  soup.faces.reserve(mesh.number_of_faces());
  for (const auto f : mesh.faces()) {
    std::vector<std::size_t> face;
    for (const auto v : CGAL::vertices_around_face(mesh.halfedge(f), mesh)) {
      face.push_back(static_cast<std::size_t>(v));
    }
    soup.faces.push_back(std::move(face));
  }
  return soup;
}

Topology topology(const SurfaceMesh& mesh) {
  ORIGAMIZER_CHECK(!mesh.has_garbage(), "the mesh has removed elements; collect garbage first");
  Topology t;
  t.vertices = mesh.number_of_vertices();
  t.edges = mesh.number_of_edges();
  t.faces = mesh.number_of_faces();

  std::vector<char> seen(mesh.number_of_halfedges(), 0);
  for (const auto h : mesh.halfedges()) {
    if (!mesh.is_border(h) || seen[static_cast<std::size_t>(h)]) continue;
    ++t.boundary_loops;
    auto g = h;
    do {
      seen[static_cast<std::size_t>(g)] = 1;
      g = mesh.next(g);
    } while (g != h);
  }

  DisjointSets joined(t.faces);
  for (const auto e : mesh.edges()) {
    if (mesh.is_border(e)) continue;
    const auto h = mesh.halfedge(e);
    joined.unite(static_cast<std::size_t>(mesh.face(h)),
                 static_cast<std::size_t>(mesh.face(mesh.opposite(h))));
  }
  std::set<std::size_t> roots;
  for (std::size_t f = 0; f < t.faces; ++f) roots.insert(joined.find(f));
  for (const auto v : mesh.vertices()) {
    ORIGAMIZER_CHECK(!mesh.is_isolated(v), "the mesh has a vertex in no face");
  }
  t.components = roots.size();

  t.euler_characteristic =
      static_cast<long>(t.vertices) - static_cast<long>(t.edges) + static_cast<long>(t.faces);
  const long twice_genus = 2 * static_cast<long>(t.components) - t.euler_characteristic -
                           static_cast<long>(t.boundary_loops);
  ORIGAMIZER_CHECK(twice_genus >= 0 && twice_genus % 2 == 0,
                   "the counts do not describe an orientable surface");
  t.genus = twice_genus / 2;
  return t;
}

double corner_angle_sum(const SurfaceMesh& mesh, SurfaceMesh::Vertex_index v) {
  const Point3& p = mesh.point(v);
  double sum = 0;
  for (const auto h : CGAL::halfedges_around_target(mesh.halfedge(v), mesh)) {
    if (mesh.is_border(h)) continue;
    const Vector3 a = mesh.point(mesh.source(h)) - p;
    const Vector3 b = mesh.point(mesh.target(mesh.next(h))) - p;
    sum += std::atan2(std::sqrt(CGAL::cross_product(a, b).squared_length()), a * b);
  }
  return sum;
}

double gaussian_curvature(const SurfaceMesh& mesh, SurfaceMesh::Vertex_index v) {
  const double full = mesh.is_border(v) ? std::numbers::pi : 2 * std::numbers::pi;
  return full - corner_angle_sum(mesh, v);
}

double signed_volume(const SurfaceMesh& mesh) {
  ORIGAMIZER_CHECK(topology(mesh).closed(), "signed_volume needs a closed mesh");
  Vector3 c(0, 0, 0);
  for (const auto v : mesh.vertices()) c = c + (mesh.point(v) - CGAL::ORIGIN);
  c = c / static_cast<double>(mesh.number_of_vertices());
  double six_volume = 0;
  for (const auto f : mesh.faces()) {
    std::vector<Vector3> corners;
    for (const auto v : CGAL::vertices_around_face(mesh.halfedge(f), mesh)) {
      corners.push_back((mesh.point(v) - CGAL::ORIGIN) - c);
    }
    for (std::size_t k = 1; k + 1 < corners.size(); ++k) {
      six_volume += corners[0] * CGAL::cross_product(corners[k], corners[k + 1]);
    }
  }
  return six_volume / 6;
}

bool orient_clean_side_out(SurfaceMesh& mesh) {
  if (!topology(mesh).closed()) return false;
  const double volume = signed_volume(mesh);
  const double extent = bounding_box_extent(mesh);
  ORIGAMIZER_CHECK(std::abs(volume) > 1e-12 * extent * extent * extent,
                   "the closed mesh encloses no volume, so its outside is undefined");
  if (volume > 0) return false;
  CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);
  return true;
}

// Pipeline -------------------------------------------------------------------

PreparedMesh prepare_mesh(PolygonSoup soup, const PrepareOptions& options) {
  PreparedMesh prepared;
  const SoupReport report = check_soup(soup);
  if (!report.ok()) {
    const bool fixable = report.ok_except_orientation() && report.orientable;
    if (fixable && options.make_orientation_consistent) {
      prepared.reoriented_faces = make_orientation_consistent(soup);
    } else {
      std::string message =
          "the input is not one connected, consistently oriented 2-manifold:\n" + report.describe();
      if (fixable) {
        message +=
            "(it is orientable: set make_orientation_consistent to reverse the disagreeing "
            "faces)\n";
      }
      throw InvariantError(message);
    }
  }
  SplitResult split = split_nonconvex_faces(soup, options.tolerances);
  prepared.face_origin = std::move(split.origin);
  prepared.split_faces = std::move(split.split);
  prepared.mesh = build_mesh(split.soup);
  prepared.flipped_to_clean_side = orient_clean_side_out(prepared.mesh);
  return prepared;
}

}  // namespace origamizer
