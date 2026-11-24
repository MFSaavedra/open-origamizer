# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Matías Saavedra

"""End-to-end tests of the M2 layout stage through the Python bindings."""

import math
from pathlib import Path

import pytest

import origamizer as og

MESHES = Path(__file__).resolve().parents[2] / "examples" / "meshes"
EXAMPLES = ["tetrahedron", "cube", "octahedron", "saddle", "torus", "cylinder", "l_prism"]
# With the rotation step's angles, nonnegative widths exist on these examples
# but not on the cube or the L prism (TVCG10 Sec. 6.4; docs/notes.tex, M2).
TRANSLATABLE = ["tetrahedron", "octahedron", "saddle", "torus", "cylinder"]


def disk_of(name):
    return og.cut_to_disk(og.load_mesh(MESHES / f"{name}.obj").mesh).disk


@pytest.mark.parametrize("name", EXAMPLES)
def test_molecule_mesh_counts(name):
    disk = disk_of(name)
    mesh = og.build_molecule_mesh(disk)
    t = disk.topology()
    boundary = sum(star.boundary for star in mesh.stars)
    assert len(mesh.stars) == t.vertices
    assert mesh.faces == t.faces
    assert mesh.digon_molecules == mesh.boundary_molecules == boundary
    assert mesh.edge_molecules == t.edges - boundary
    assert mesh.polygons == t.faces + boundary
    kinds = [m.kind for m in mesh.molecules]
    assert kinds.count(og.MoleculeKind.EDGE) == mesh.edge_molecules
    # Gauss-Bonnet on the disk: the defects sum to 2 pi plus pi per boundary vertex.
    assert math.isclose(sum(s.defect() for s in mesh.stars), 2 * math.pi + math.pi * boundary)


@pytest.mark.parametrize("name", EXAMPLES)
def test_mapping_and_reconstruction(name):
    disk = disk_of(name)
    mesh = og.build_molecule_mesh(disk)
    rotation = og.solve_rotation(mesh)
    assert rotation.feasible, rotation.message
    scale = og.mean_molecule_length(mesh)
    for objective in og.WidthObjective.MINIMUM_WIDTH, og.WidthObjective.PAPER_PENALTY:
        translation = og.solve_translation(mesh, rotation.theta, objective=objective)
        assert translation.feasible == (name in TRANSLATABLE), translation.message
    # Objective B returns widths that solve T17 even when some are negative.
    widths = og.solve_translation(mesh, rotation.theta,
                                  objective=og.WidthObjective.PAPER_PENALTY).w
    layout = og.reconstruct_layout(disk, mesh, rotation.theta, widths)
    assert layout.closure_error < 1e-9 * scale
    assert len(layout.corners) == mesh.polygons
    assert all(len(c) == len(v) for c, v in zip(layout.corners, layout.vertices))


def test_infeasible_translation_reports_the_violation():
    mesh = og.build_molecule_mesh(disk_of("cube"))
    theta = og.solve_rotation(mesh).theta
    result = og.solve_translation(mesh, theta, objective=og.WidthObjective.PAPER_PENALTY)
    assert not result.feasible
    assert "T5 and T7 cannot hold" in result.message
    assert math.isclose(result.violation, sum(max(0.0, -w) for w in result.w), rel_tol=1e-9)


def test_translation_options():
    mesh = og.build_molecule_mesh(disk_of("octahedron"))
    theta = og.solve_rotation(mesh).theta
    wide = og.solve_translation(mesh, theta, minimum_width=0.2)
    assert wide.feasible and min(wide.w) > 0.2 - 1e-9
    closest = og.solve_translation(mesh, theta, objective=og.WidthObjective.PAPER_PENALTY,
                                   reference=wide.w)
    assert closest.feasible
    assert max(abs(a - b) for a, b in zip(closest.w, wide.w)) < 1e-9


def test_rotation_without_digons_reports_why():
    mesh = og.build_molecule_mesh(disk_of("cube"), boundary_digons=False)
    rotation = og.solve_rotation(mesh)
    assert not rotation.feasible
    assert "need at least three" in rotation.message


def test_reconstruction_rejects_inconsistent_widths():
    disk = disk_of("octahedron")
    mesh = og.build_molecule_mesh(disk)
    theta = og.solve_rotation(mesh).theta
    w = og.solve_translation(mesh, theta).w
    w[0] += 0.1
    with pytest.raises(og.InvariantError, match="does not close"):
        og.reconstruct_layout(disk, mesh, theta, w)


def test_remapping_makes_the_cube_translatable():
    mesh = og.build_molecule_mesh(disk_of("cube"))
    plain = og.solve_rotation(mesh)
    assert not og.solve_translation(mesh, plain.theta).feasible
    remapped = og.solve_rotation(mesh, reference=og.remap_angles(mesh))
    assert remapped.feasible
    assert og.solve_translation(mesh, remapped.theta).feasible


def test_the_l_prism_needs_the_whole_cut_tree():
    # With the minimal cut no nonnegative widths exist, even after remapping;
    # with every vertex on the cut they do (docs/notes.tex, M2).
    mesh = og.load_mesh(MESHES / "l_prism.obj").mesh
    for strategy, expected in (og.CutStrategy.MINIMAL, False), (og.CutStrategy.TREE, True):
        molecules = og.build_molecule_mesh(og.cut_to_disk(mesh, strategy=strategy).disk)
        rotation = og.solve_rotation(molecules, reference=og.remap_angles(molecules))
        assert rotation.feasible
        assert og.solve_translation(molecules, rotation.theta).feasible == expected
