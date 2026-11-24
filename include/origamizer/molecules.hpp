// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

/// \file
/// The molecule mesh of TVCG10 on a surface cut into a disk (milestone M2):
/// one edge-tucking molecule per interior edge, one boundary molecule per
/// boundary vertex, optionally the boundary digons of TVCG10 Sec. 5.1, and
/// the star of every vertex, which lists its molecules and sector angles
/// counterclockwise (TVCG10 Sec. 3 and 4.1; equation tags refer to
/// docs/tvcg2010-notes; design notes in docs/notes.tex, M2).

#pragma once

#include <cstddef>
#include <vector>

#include "origamizer/mesh.hpp"

namespace origamizer {

/// A molecule and its angle variable theta, seen from vertex i.
struct Molecule {
  enum class Kind {
    /// ETM(i, j) between the two face copies of an interior edge (TVCG10 Sec. 3.1).
    kEdge,
    /// ETM between a boundary face and the digon added along its boundary
    /// edge (i, j) (TVCG10 Sec. 5.1).
    kDigon,
    /// ETM(i, o): the outside of the paper at boundary vertex i (TVCG10
    /// Sec. 4.1; our reading in docs/tvcg2010-notes, Sec. 4).
    kBoundary,
  };
  Kind kind = Kind::kEdge;
  std::size_t i = 0;  ///< vertex the molecule is oriented from
  std::size_t j = 0;  ///< other vertex; equals i for a boundary molecule
  double length = 0;  ///< l(i, j); 0 for a boundary molecule
};

/// One slot j_n of a vertex star.
struct StarSlot {
  std::size_t molecule = 0;  ///< index into MoleculeMesh::molecules
  int sign = 1;              ///< theta(v, j_n) = sign * theta of the molecule
  /// alpha(v, j_n): the corner angle of the polygon between slot n-1 and
  /// slot n (TVCG10 Eq. 1); 0 for a digon corner.
  double sector = 0;
  /// The polygon between slot n-1 and slot n, whose corner `sector` is: a
  /// face index of the disk, or for a digon, MoleculeMesh::faces plus the
  /// digon's index among the digon molecules.
  std::size_t polygon = 0;
};

/// The molecules around a vertex in counterclockwise order, seen from the
/// clean side. A boundary vertex has exactly one boundary-molecule slot.
struct VertexStar {
  std::vector<StarSlot> slots;
  bool boundary = false;

  /// 2 pi minus the sum of the sector angles: the right side of TVCG10 Eq. 1
  /// (T1), which the molecules around the vertex must absorb.
  double defect() const;
};

struct MoleculeOptions {
  /// Add the boundary molecules of TVCG10 Sec. 5.1: along every boundary
  /// edge, a zero-area digon glued to its face by an extra edge-tucking
  /// molecule. Every boundary star gains two slots. Without them, a boundary
  /// vertex with fewer than three faces has fewer than three slots, and T1
  /// with T8 cannot hold there (docs/notes.tex, M2).
  bool boundary_digons = true;
};

struct MoleculeMesh {
  /// Edge molecules first, then digon molecules, then boundary molecules.
  std::vector<Molecule> molecules;
  std::vector<VertexStar> stars;  ///< one per vertex of the disk
  std::size_t edge_molecules = 0;
  std::size_t digon_molecules = 0;
  std::size_t faces = 0;  ///< faces of the disk: the polygons before the digons

  std::size_t boundary_molecules() const {
    return molecules.size() - edge_molecules - digon_molecules;
  }
  /// Faces, then one digon per digon molecule.
  std::size_t polygons() const { return faces + digon_molecules; }
};

/// Builds the molecule mesh of a topological disk with consistently oriented
/// faces (TVCG10 Sec. 3, 4.1, and 5.1). Edge and digon molecules are oriented
/// from their smaller vertex index to the larger one, in edge order; boundary
/// molecules follow, one per boundary vertex in vertex order. Throws
/// InvariantError unless `disk` is a disk.
MoleculeMesh build_molecule_mesh(const SurfaceMesh& disk, const MoleculeOptions& options = {});

}  // namespace origamizer
