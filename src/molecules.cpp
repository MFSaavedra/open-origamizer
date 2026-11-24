// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

#include "origamizer/molecules.hpp"

#include <cmath>
#include <limits>
#include <numbers>
#include <string>

#include "origamizer/check.hpp"

namespace origamizer {
namespace {

using Halfedge = SurfaceMesh::Halfedge_index;
using Vertex = SurfaceMesh::Vertex_index;

constexpr std::size_t kNone = std::numeric_limits<std::size_t>::max();

template <class Index>
std::size_t to_size(Index i) {
  return static_cast<std::size_t>(i);
}

/// Corner angle of face(h) at target(h), for a non-border halfedge h.
double corner_angle(const SurfaceMesh& mesh, Halfedge h) {
  const Point3& p = mesh.point(mesh.target(h));
  const Vector3 a = mesh.point(mesh.source(h)) - p;
  const Vector3 b = mesh.point(mesh.target(mesh.next(h))) - p;
  return std::atan2(std::sqrt(CGAL::cross_product(a, b).squared_length()), a * b);
}

/// The halfedges into v whose faces surround v, in counterclockwise order.
///
/// With counterclockwise faces, the corner of face(h) at v lies between
/// edge(next(h)) on its clockwise side and edge(h) on its counterclockwise
/// side, and the next face counterclockwise is the face of prev(opposite(h)).
/// (CGAL's halfedges_around_target runs the other way.) At a boundary vertex
/// the sequence starts at the face whose clockwise side is a boundary edge, so
/// the last halfedge lies on the other boundary edge.
std::vector<Halfedge> counterclockwise_corners(const SurfaceMesh& mesh, Vertex v) {
  Halfedge start = SurfaceMesh::null_halfedge();
  for (const Halfedge h : CGAL::halfedges_around_target(v, mesh)) {
    if (mesh.is_border(h)) continue;
    if (start == SurfaceMesh::null_halfedge()) start = h;
    if (mesh.is_border(mesh.opposite(mesh.next(h)))) {  // clockwise side is on the boundary
      start = h;
      break;
    }
  }
  ORIGAMIZER_CHECK(start != SurfaceMesh::null_halfedge(), "a vertex has no face");
  std::vector<Halfedge> corners;
  Halfedge h = start;
  do {
    corners.push_back(h);
    const Halfedge o = mesh.opposite(h);
    if (mesh.is_border(o)) break;  // reached the counterclockwise boundary edge
    h = mesh.prev(o);
  } while (h != start);
  return corners;
}

}  // namespace

double VertexStar::defect() const {
  double sum = 0;
  for (const StarSlot& slot : slots) sum += slot.sector;
  return 2 * std::numbers::pi - sum;
}

MoleculeMesh build_molecule_mesh(const SurfaceMesh& disk, const MoleculeOptions& options) {
  const Topology t = topology(disk);
  ORIGAMIZER_CHECK(t.is_disk(),
                   "the molecule mesh needs a topological disk; cut the surface first");

  MoleculeMesh mesh;
  // Edge molecules (interior edges), then digon molecules (boundary edges).
  std::vector<std::size_t> edge_molecule(disk.number_of_edges(), kNone);
  for (const bool border : {false, true}) {
    if (border && !options.boundary_digons) break;
    for (const auto e : disk.edges()) {
      if (disk.is_border(e) != border) continue;
      const Halfedge h = disk.halfedge(e);
      const std::size_t a = to_size(disk.source(h));
      const std::size_t b = to_size(disk.target(h));
      Molecule molecule;
      molecule.kind = border ? Molecule::Kind::kDigon : Molecule::Kind::kEdge;
      molecule.i = std::min(a, b);
      molecule.j = std::max(a, b);
      molecule.length =
          std::sqrt(CGAL::squared_distance(disk.point(disk.source(h)), disk.point(disk.target(h))));
      edge_molecule[to_size(e)] = mesh.molecules.size();
      mesh.molecules.push_back(molecule);
    }
    (border ? mesh.digon_molecules : mesh.edge_molecules) =
        mesh.molecules.size() - mesh.edge_molecules;
  }

  mesh.stars.resize(t.vertices);
  std::vector<int> appearances(mesh.molecules.size(), 0);
  const auto edge_slot = [&](Halfedge h, std::size_t v, double sector) {
    StarSlot slot;
    slot.molecule = edge_molecule[to_size(disk.edge(h))];
    slot.sign = mesh.molecules[slot.molecule].i == v ? 1 : -1;
    slot.sector = sector;
    ++appearances[slot.molecule];
    return slot;
  };
  for (const Vertex v : disk.vertices()) {
    VertexStar& star = mesh.stars[to_size(v)];
    star.boundary = disk.is_border(v);
    const std::vector<Halfedge> corners = counterclockwise_corners(disk, v);
    for (const Halfedge h : corners) {
      if (!disk.is_border(disk.edge(h))) {
        star.slots.push_back(edge_slot(h, to_size(v), corner_angle(disk, h)));
        continue;
      }
      // The counterclockwise boundary edge. Going on counterclockwise: the
      // digon molecule of that edge, the outside (slot o), and the digon
      // molecule of the clockwise boundary edge, across two zero-angle digon
      // corners; without digons, slot o alone.
      StarSlot outside;
      outside.molecule = mesh.molecules.size();
      Molecule molecule;
      molecule.kind = Molecule::Kind::kBoundary;
      molecule.i = molecule.j = to_size(v);
      mesh.molecules.push_back(molecule);
      appearances.push_back(1);
      if (options.boundary_digons) {
        star.slots.push_back(edge_slot(h, to_size(v), corner_angle(disk, h)));
        outside.sector = 0;
        star.slots.push_back(outside);
        star.slots.push_back(edge_slot(disk.next(corners.front()), to_size(v), 0));
      } else {
        outside.sector = corner_angle(disk, h);
        star.slots.push_back(outside);
      }
    }
    ORIGAMIZER_CHECK(!star.slots.empty(), "a vertex star is empty");
  }

  // Every edge or digon molecule sits in the stars of both its ends with
  // opposite signs; every boundary vertex has exactly one boundary molecule.
  for (std::size_t n = 0; n < mesh.molecules.size(); ++n) {
    const bool shared = mesh.molecules[n].kind != Molecule::Kind::kBoundary;
    ORIGAMIZER_CHECK(appearances[n] == (shared ? 2 : 1),
                     "molecule " + std::to_string(n) + " appears in the wrong number of stars");
  }
  for (std::size_t v = 0; v < mesh.stars.size(); ++v) {
    int boundary_slots = 0;
    for (const StarSlot& slot : mesh.stars[v].slots) {
      boundary_slots += mesh.molecules[slot.molecule].kind == Molecule::Kind::kBoundary;
    }
    ORIGAMIZER_CHECK(boundary_slots == (mesh.stars[v].boundary ? 1 : 0),
                     "vertex " + std::to_string(v) + " has the wrong number of boundary molecules");
  }
  return mesh;
}

}  // namespace origamizer
