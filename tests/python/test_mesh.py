# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Matías Saavedra

"""End-to-end tests of the M1 mesh stage through the Python bindings."""

import math
from pathlib import Path

import pytest

import origamizer as og

MESHES = Path(__file__).resolve().parents[2] / "examples" / "meshes"
EXAMPLES = ["tetrahedron", "cube", "octahedron", "saddle", "torus", "cylinder", "l_prism"]

# (Euler characteristic, boundary loops, genus) of each example.
TOPOLOGY = {
    "tetrahedron": (2, 0, 0),
    "cube": (2, 0, 0),
    "octahedron": (2, 0, 0),
    "saddle": (1, 1, 0),
    "torus": (0, 0, 1),
    "cylinder": (0, 2, 0),
    "l_prism": (2, 0, 0),
}


def load(name):
    return og.load_mesh(MESHES / f"{name}.obj")


@pytest.mark.parametrize("name", EXAMPLES)
def test_example_topology_and_gauss_bonnet(name):
    mesh = load(name).mesh
    t = mesh.topology()
    assert (t.euler_characteristic, t.boundary_loops, t.genus) == TOPOLOGY[name]
    assert math.isclose(sum(mesh.gaussian_curvature()), 2 * math.pi * t.euler_characteristic,
                        abs_tol=1e-9)


@pytest.mark.parametrize("name", EXAMPLES)
def test_example_cuts_into_a_disk(name):
    mesh = load(name).mesh
    cut = og.cut_to_disk(mesh)
    assert cut.disk.topology().is_disk
    assert len(cut.disk.faces) == len(mesh.faces)
    points = mesh.points
    for v, original in enumerate(cut.original_vertex):
        assert cut.disk.points[v] == points[original]
    assert math.isclose(sum(cut.disk.gaussian_curvature()), 2 * math.pi, abs_tol=1e-9)


def test_cube_cut_and_boundary():
    cut = og.cut_to_disk(load("cube").mesh)
    assert cut.cut_edges == [(0, 1), (1, 3), (3, 7)]
    loops = cut.disk.boundary_loops()
    assert len(loops) == 1 and len(loops[0]) == 6


def test_explicit_cut_edges():
    mesh = load("tetrahedron").mesh
    cut = og.cut_to_disk(mesh, cut_edges=[(0, 2), (2, 3)])
    assert cut.cut_edges == [(0, 2), (2, 3)]
    with pytest.raises(og.InvariantError, match="would be glued back"):
        og.cut_along(mesh, [(1, 3)])


def test_l_prism_caps_are_split():
    prepared = load("l_prism")
    assert prepared.split_faces == [0, 1]
    assert len(prepared.mesh.faces) == 10
    assert prepared.face_origin[:4] == [0, 0, 1, 1]


def test_soup_round_trip_and_checks():
    soup = og.PolygonSoup([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 1]], [[0, 1, 2], [0, 1, 3]])
    report = og.check_soup(soup)
    assert not report.consistently_oriented and report.orientable
    assert report.ok_except_orientation()
    with pytest.raises(og.InvariantError, match="make_orientation_consistent"):
        og.prepare_mesh(soup)
    prepared = og.prepare_mesh(soup, make_orientation_consistent=True)
    assert prepared.reoriented_faces == [1]


def test_parse_errors_name_the_line():
    with pytest.raises(og.InvariantError, match=r"inline\.obj:2: unsupported OBJ element 'l'"):
        og.parse_obj("v 0 0 0\nl 1 2\n", "inline.obj")


def test_missing_file_is_a_runtime_error_not_an_invariant_error():
    with pytest.raises(RuntimeError) as info:
        og.read_obj(MESHES / "does-not-exist.obj")
    assert not isinstance(info.value, og.InvariantError)


def test_nonplanar_face_is_rejected(tmp_path):
    path = tmp_path / "warped.obj"
    path.write_text("v 0 0 0\nv 1 0 0\nv 1 1 0.3\nv 0 1 0\nf 1 2 3 4\n")
    with pytest.raises(og.InvariantError, match="not planar"):
        og.load_mesh(path)
