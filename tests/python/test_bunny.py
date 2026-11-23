# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Matías Saavedra

"""The low-poly Stanford bunny, generated locally by examples/meshes/fetch_bunny.py."""

import math
from pathlib import Path

import pytest

import origamizer as og

BUNNY = Path(__file__).resolve().parents[2] / "examples" / "meshes" / "bunny.obj"

pytestmark = pytest.mark.skipif(
    not BUNNY.exists(), reason="run examples/meshes/fetch_bunny.py to generate the bunny")


def test_bunny_prepares_and_cuts_into_a_disk():
    mesh = og.load_mesh(BUNNY).mesh
    t = mesh.topology()
    assert t.genus == 0
    assert t.faces <= 374
    assert math.isclose(sum(mesh.gaussian_curvature()), 2 * math.pi * t.euler_characteristic,
                        abs_tol=1e-9)
    cut = og.cut_to_disk(mesh)
    assert cut.disk.topology().is_disk
    # The holes of the scan are joined by paths: no cut vertex is a dangling tip.
    degree = {}
    for a, b in cut.cut_edges:
        degree[a] = degree.get(a, 0) + 1
        degree[b] = degree.get(b, 0) + 1
    boundary = {v for loop in mesh.boundary_loops() for v in loop}
    assert all(d >= 2 or v in boundary for v, d in degree.items())
