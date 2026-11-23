// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

#include "origamizer/cut.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <string>
#include <utility>

#include "detail/disjoint_sets.hpp"
#include "origamizer/check.hpp"

namespace origamizer {
namespace {

using Vertex = SurfaceMesh::Vertex_index;
using Halfedge = SurfaceMesh::Halfedge_index;
using Edge = SurfaceMesh::Edge_index;
using Face = SurfaceMesh::Face_index;

constexpr std::size_t kNone = std::numeric_limits<std::size_t>::max();

template <class Index>
std::size_t to_size(Index i) {
  return static_cast<std::size_t>(i);
}

Vertex vertex_at(std::size_t i) { return Vertex(static_cast<SurfaceMesh::size_type>(i)); }

Face face_at(std::size_t i) { return Face(static_cast<SurfaceMesh::size_type>(i)); }

IndexPair edge_vertices(const SurfaceMesh& mesh, Edge e) {
  const Halfedge h = mesh.halfedge(e);
  const std::size_t a = to_size(mesh.source(h));
  const std::size_t b = to_size(mesh.target(h));
  return {std::min(a, b), std::max(a, b)};
}

std::string format_pair(const IndexPair& p) {
  return "(" + std::to_string(p.first) + ", " + std::to_string(p.second) + ")";
}

/// Shortest edge-path distances from `source` (Dijkstra on edge lengths) and,
/// for each vertex, the halfedge that reaches it on a shortest path.
struct ShortestPaths {
  std::vector<double> distance;
  std::vector<Halfedge> via;
};

ShortestPaths shortest_paths(const SurfaceMesh& mesh, Vertex source) {
  const std::size_t n = mesh.number_of_vertices();
  ShortestPaths paths{std::vector<double>(n, std::numeric_limits<double>::infinity()),
                      std::vector<Halfedge>(n, SurfaceMesh::null_halfedge())};
  using Item = std::pair<double, std::size_t>;
  std::priority_queue<Item, std::vector<Item>, std::greater<>> pending;
  paths.distance[to_size(source)] = 0;
  pending.push({0.0, to_size(source)});
  while (!pending.empty()) {
    const auto [d, u] = pending.top();
    pending.pop();
    if (d > paths.distance[u]) continue;
    for (const Halfedge h : CGAL::halfedges_around_target(vertex_at(u), mesh)) {
      const Vertex w = mesh.source(h);  // h runs from w to u
      const double next =
          d + std::sqrt(CGAL::squared_distance(mesh.point(w), mesh.point(vertex_at(u))));
      if (next < paths.distance[to_size(w)]) {
        paths.distance[to_size(w)] = next;
        paths.via[to_size(w)] = mesh.opposite(h);  // from u to w
        pending.push({next, to_size(w)});
      }
    }
  }
  return paths;
}

/// One shortest edge path between the two most curved vertices of a closed
/// genus-0 surface (see automatic_cut).
std::vector<IndexPair> seed_path(const SurfaceMesh& mesh) {
  const std::size_t n = mesh.number_of_vertices();
  ORIGAMIZER_CHECK(n >= 2, "a closed surface needs at least two vertices");
  std::vector<double> curvature(n);
  for (const Vertex v : mesh.vertices()) curvature[to_size(v)] = gaussian_curvature(mesh, v);
  constexpr double kTie = 1e-9;  // curvatures closer than this count as equal
  const auto more_curved = [&curvature](std::size_t a, std::size_t b) {
    return curvature[a] > curvature[b] + kTie;
  };

  std::size_t a = 0;
  for (std::size_t v = 1; v < n; ++v) {
    if (more_curved(v, a)) a = v;
  }
  const ShortestPaths paths = shortest_paths(mesh, vertex_at(a));
  std::size_t b = kNone;
  for (std::size_t v = 0; v < n; ++v) {
    if (v == a) continue;
    if (b == kNone || more_curved(v, b) ||
        (!more_curved(b, v) && paths.distance[v] > paths.distance[b])) {
      b = v;
    }
  }
  ORIGAMIZER_CHECK(std::isfinite(paths.distance[b]), "the surface is not connected");

  std::vector<IndexPair> path;
  for (Vertex v = vertex_at(b); v != vertex_at(a);) {
    const Halfedge h = paths.via[to_size(v)];
    path.push_back(edge_vertices(mesh, mesh.edge(h)));
    v = mesh.source(h);
  }

  // A one-edge slit leaves both ends unsplit, so its two sides would be glued
  // back; extend it from b to the best neighbor by the same criteria.
  if (path.size() == 1) {
    std::size_t c = kNone;
    for (const Halfedge h : CGAL::halfedges_around_target(vertex_at(b), mesh)) {
      const std::size_t w = to_size(mesh.source(h));
      if (w == a) continue;
      if (c == kNone || more_curved(w, c) ||
          (!more_curved(c, w) && (paths.distance[w] > paths.distance[c] ||
                                  (paths.distance[w] == paths.distance[c] && w < c)))) {
        c = w;
      }
    }
    ORIGAMIZER_CHECK(c != kNone, "a vertex of a closed surface has only one neighbor");
    path.push_back({std::min(b, c), std::max(b, c)});
  }
  return path;
}

}  // namespace

std::vector<IndexPair> automatic_cut(const SurfaceMesh& mesh, std::size_t root_face) {
  const Topology t = topology(mesh);
  ORIGAMIZER_CHECK(t.components == 1, "only a connected surface can be cut into one disk");
  ORIGAMIZER_CHECK(root_face < t.faces, "root face " + std::to_string(root_face) +
                                            " does not exist; the mesh has " +
                                            std::to_string(t.faces) + " faces");

  // Breadth-first flooding (TVCG10 Sec. 8.1): crossed edges stay glued.
  std::vector<char> glued(mesh.number_of_edges(), 0);
  std::vector<char> reached(t.faces, 0);
  std::queue<Face> pending;
  reached[root_face] = 1;
  pending.push(face_at(root_face));
  while (!pending.empty()) {
    const Face f = pending.front();
    pending.pop();
    for (const Halfedge h : CGAL::halfedges_around_face(mesh.halfedge(f), mesh)) {
      const Halfedge o = mesh.opposite(h);
      if (mesh.is_border(o)) continue;
      const Face g = mesh.face(o);
      if (reached[to_size(g)]) continue;
      reached[to_size(g)] = 1;
      glued[to_size(mesh.edge(h))] = 1;
      pending.push(g);
    }
  }

  // Cut graph: the interior edges that were not crossed, plus the boundary,
  // which is never cut but anchors the pruning (boundary vertices keep degree
  // at least 2).
  std::vector<char> in_graph(mesh.number_of_edges(), 0);
  std::vector<std::size_t> degree(t.vertices, 0);
  for (const Edge e : mesh.edges()) {
    if (glued[to_size(e)]) continue;
    in_graph[to_size(e)] = 1;
    const auto [a, b] = edge_vertices(mesh, e);
    ++degree[a];
    ++degree[b];
  }

  // Pruning: glue back the single cut edge at the tip of a slit.
  std::queue<Vertex> tips;
  for (const Vertex v : mesh.vertices()) {
    if (degree[to_size(v)] == 1) tips.push(v);
  }
  while (!tips.empty()) {
    const Vertex v = tips.front();
    tips.pop();
    if (degree[to_size(v)] != 1) continue;
    for (const Halfedge h : CGAL::halfedges_around_target(v, mesh)) {
      const Edge e = mesh.edge(h);
      if (!in_graph[to_size(e)]) continue;
      ORIGAMIZER_CHECK(!mesh.is_border(e), "pruning reached a boundary edge");
      in_graph[to_size(e)] = 0;
      --degree[to_size(v)];
      const Vertex w = mesh.source(h);
      if (--degree[to_size(w)] == 1) tips.push(w);
      break;
    }
  }

  std::vector<IndexPair> cut;
  for (const Edge e : mesh.edges()) {
    if (in_graph[to_size(e)] && !mesh.is_border(e)) cut.push_back(edge_vertices(mesh, e));
  }
  if (cut.empty() && t.closed()) {
    ORIGAMIZER_CHECK(t.genus == 0, "pruning removed every cut of a surface with handles");
    cut = seed_path(mesh);
  }
  std::sort(cut.begin(), cut.end());
  return cut;
}

CutResult cut_along(const SurfaceMesh& mesh, const std::vector<IndexPair>& edges) {
  ORIGAMIZER_CHECK(!mesh.has_garbage(), "the mesh has removed elements; collect garbage first");
  const std::size_t n = mesh.number_of_vertices();
  std::vector<char> cut(mesh.number_of_edges(), 0);
  for (const auto& [a, b] : edges) {
    const std::string name = format_pair({a, b});
    ORIGAMIZER_CHECK(a < n && b < n && a != b,
                     "cut edge " + name + " does not name two vertices of the mesh");
    const auto [h, found] = CGAL::halfedge(vertex_at(a), vertex_at(b), mesh);
    ORIGAMIZER_CHECK(found, name + " is not an edge of the mesh");
    ORIGAMIZER_CHECK(!mesh.is_border(mesh.edge(h)),
                     name + " is a boundary edge; only interior edges can be cut");
    cut[to_size(mesh.edge(h))] = 1;
  }

  CutResult result;
  for (const Edge e : mesh.edges()) {
    if (cut[to_size(e)]) result.cut_edges.push_back(edge_vertices(mesh, e));
  }
  std::sort(result.cut_edges.begin(), result.cut_edges.end());

  // A corner is a non-border halfedge h: the corner of face(h) at target(h).
  // Corners on both sides of a glued edge belong to the same copy of a vertex.
  detail::DisjointSets corners(mesh.number_of_halfedges());
  for (const Edge e : mesh.edges()) {
    if (cut[to_size(e)] || mesh.is_border(e)) continue;
    const Halfedge h = mesh.halfedge(e);
    const Halfedge o = mesh.opposite(h);
    corners.unite(to_size(h), to_size(mesh.prev(o)));  // corners at target(h)
    corners.unite(to_size(o), to_size(mesh.prev(h)));  // corners at source(h)
  }

  // A cut edge whose ends are both unsplit would be glued back together.
  for (const Edge e : mesh.edges()) {
    if (!cut[to_size(e)]) continue;
    const Halfedge h = mesh.halfedge(e);
    const Halfedge o = mesh.opposite(h);
    const bool target_unsplit = corners.find(to_size(h)) == corners.find(to_size(mesh.prev(o)));
    const bool source_unsplit = corners.find(to_size(o)) == corners.find(to_size(mesh.prev(h)));
    ORIGAMIZER_CHECK(!(target_unsplit && source_unsplit),
                     "cut edge " + format_pair(edge_vertices(mesh, e)) +
                         " would be glued back: neither of its ends is split by another cut "
                         "edge or by the boundary");
  }

  // Number the copies vertex by vertex, so copies of a vertex are adjacent.
  PolygonSoup soup;
  std::vector<std::size_t> copy(mesh.number_of_halfedges(), kNone);
  for (const Vertex v : mesh.vertices()) {
    for (const Halfedge h : CGAL::halfedges_around_target(v, mesh)) {
      if (mesh.is_border(h)) continue;
      const std::size_t root = corners.find(to_size(h));
      if (copy[root] != kNone) continue;
      copy[root] = soup.points.size();
      soup.points.push_back(mesh.point(v));
      result.original_vertex.push_back(to_size(v));
    }
  }
  for (const Face f : mesh.faces()) {
    std::vector<std::size_t> face;
    for (const Halfedge h : CGAL::halfedges_around_face(mesh.halfedge(f), mesh)) {
      face.push_back(copy[corners.find(to_size(h))]);
    }
    soup.faces.push_back(std::move(face));
  }

  try {
    result.disk = build_mesh(soup, /*allow_coincident_vertices=*/true);
  } catch (const InvariantError& error) {
    throw InvariantError("cutting along the given edges does not leave one surface:\n" +
                         std::string(error.what()));
  }
  const Topology t = topology(result.disk);
  ORIGAMIZER_CHECK(t.is_disk(), "cutting along " + std::to_string(result.cut_edges.size()) +
                                    " edges leaves a surface with Euler characteristic " +
                                    std::to_string(t.euler_characteristic) + " and " +
                                    std::to_string(t.boundary_loops) +
                                    " boundary loops, not a topological disk");
  return result;
}

CutResult cut_to_disk(const SurfaceMesh& mesh, const CutOptions& options) {
  return cut_along(
      mesh, options.cut_edges.empty() ? automatic_cut(mesh, options.root_face) : options.cut_edges);
}

}  // namespace origamizer
