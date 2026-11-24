// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

/// \file
/// Cutting a surface into a topological disk (milestone M1).
///
/// The molecule mesh of TVCG10 lives on a disk: a closed surface, or one with
/// handles or several boundary loops, is first cut along mesh edges into a
/// polygonal schema (TVCG10 Sec. 2.2.1 and 8.1; SoCG17 Sec. 2). The two
/// copies of each cut edge end up on the boundary of the paper.

#pragma once

#include <cstddef>
#include <vector>

#include "origamizer/mesh.hpp"

namespace origamizer {

/// How much of the flooding's cut graph the automatic cut keeps.
enum class CutStrategy {
  /// The essential part only, the fewest cut edges: 2g loops for genus g,
  /// paths between boundary loops, or a seed path on a closed genus-0
  /// surface (see automatic_cut).
  kMinimal,
  /// The whole cut tree: every interior edge the flooding does not cross, as
  /// in TVCG10 Sec. 8.1 before its simplification. Every vertex then lies on
  /// the cut, which gives the mapping the most freedom at the cost of the
  /// most seams (docs/notes.tex, M2).
  kTree,
};

struct CutOptions {
  /// Face where the breadth-first flooding starts. TVCG10 Sec. 8.1 lets the
  /// user choose it; it ends up near the middle of the paper.
  std::size_t root_face = 0;
  CutStrategy strategy = CutStrategy::kMinimal;
  /// Edges to cut, as pairs of vertex indices. When non-empty they replace
  /// the automatic choice and must cut the surface into a disk.
  std::vector<IndexPair> cut_edges;
};

struct CutResult {
  /// The cut surface, a topological disk. Face i is face i of the input.
  SurfaceMesh disk;
  /// For each vertex of `disk`, the input vertex it is a copy of.
  std::vector<std::size_t> original_vertex;
  /// The input edges that were cut, as sorted vertex pairs, in sorted order.
  std::vector<IndexPair> cut_edges;
};

/// Chooses interior edges whose removal leaves a topological disk.
///
/// Floods the faces breadth first from `root_face` (TVCG10 Sec. 8.1); the
/// edges the flooding crosses form a spanning tree of the dual graph and stay
/// glued, so cutting all other interior edges leaves a disk. Then, as long as
/// some vertex off the original boundary touches exactly one cut edge, that
/// edge is glued back; this closes a slit from its tip and keeps a disk. What
/// remains are 2g loops for a surface of genus g, plus paths joining its
/// boundary loops. A closed genus-0 surface prunes away completely, so it is
/// cut along one shortest edge path between its two most curved vertices
/// (largest angle defect, then largest distance, then smallest index). A path
/// of one edge is extended by one more edge chosen the same way, because a
/// one-edge slit splits neither end and would be glued back together.
/// TVCG10 Sec. 10 notes that a boundary through sharp points tends to give
/// more efficient crease patterns; finer control of the cut is future work.
/// With CutStrategy::kTree, the pruning and the seed path are skipped.
std::vector<IndexPair> automatic_cut(const SurfaceMesh& mesh, std::size_t root_face = 0,
                                     CutStrategy strategy = CutStrategy::kMinimal);

/// Cuts the mesh along `edges`: every vertex gets one copy per fan of faces
/// between cut edges. Throws InvariantError if a pair is not an interior edge
/// of the mesh, if a cut edge would be glued back (neither end split), or if
/// the result is not a topological disk.
CutResult cut_along(const SurfaceMesh& mesh, const std::vector<IndexPair>& edges);

/// automatic_cut, or the edges in `options`, followed by cut_along.
CutResult cut_to_disk(const SurfaceMesh& mesh, const CutOptions& options = {});

}  // namespace origamizer
