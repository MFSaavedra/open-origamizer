# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Matías Saavedra

"""Open Origamizer: fold one uncut sheet into a polyhedral surface by tucking.

The algorithms live in the C++ extension module ``origamizer._core``. This
package adds the file writers, the fabrication helpers, and the command line.
"""

from os import PathLike

from ._core import (
    CutResult,
    CutStrategy,
    FaceShape,
    InvariantError,
    Layout,
    Mesh,
    Molecule,
    MoleculeKind,
    MoleculeMesh,
    PolygonSoup,
    PreparedMesh,
    RotationResult,
    RotationViolation,
    RotationViolationKind,
    SoupReport,
    StarSlot,
    Topology,
    TranslationResult,
    VertexStar,
    WidthObjective,
    automatic_cut,
    build_mesh,
    build_molecule_mesh,
    cgal_version,
    check_soup,
    cut_along,
    cut_to_disk,
    exterior_angles,
    face_shape,
    mean_molecule_length,
    parse_obj,
    prepare_mesh,
    read_obj,
    reconstruct_layout,
    remap_angles,
    simplify_triangle_mesh,
    solve_rotation,
    solve_translation,
    to_soup,
)


def load_mesh(path: str | PathLike, **options) -> PreparedMesh:
    """Read an OBJ file and prepare it; keyword options go to prepare_mesh."""
    return prepare_mesh(read_obj(path), **options)


__all__ = [
    "CutResult",
    "CutStrategy",
    "FaceShape",
    "InvariantError",
    "Layout",
    "Mesh",
    "Molecule",
    "MoleculeKind",
    "MoleculeMesh",
    "PolygonSoup",
    "PreparedMesh",
    "RotationResult",
    "RotationViolation",
    "RotationViolationKind",
    "SoupReport",
    "StarSlot",
    "Topology",
    "TranslationResult",
    "VertexStar",
    "WidthObjective",
    "automatic_cut",
    "build_mesh",
    "build_molecule_mesh",
    "cgal_version",
    "check_soup",
    "cut_along",
    "cut_to_disk",
    "exterior_angles",
    "face_shape",
    "load_mesh",
    "mean_molecule_length",
    "parse_obj",
    "prepare_mesh",
    "read_obj",
    "reconstruct_layout",
    "remap_angles",
    "simplify_triangle_mesh",
    "solve_rotation",
    "solve_translation",
    "to_soup",
]
