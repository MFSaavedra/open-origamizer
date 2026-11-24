# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Matías Saavedra

"""The sample models of the original Origamizer, when copied locally.

See examples/meshes/origamizer/README.md: the files are not committed, so these
tests skip when they are absent (always in CI). Results and their discussion
are in docs/notes.tex, M2.
"""

import math
from pathlib import Path

import pytest

import origamizer as og

SAMPLES = Path(__file__).resolve().parents[2] / "examples" / "meshes" / "origamizer"
# (faces after triangulating nonplanar quads, quads triangulated)
EXPECTED = {
    "bell8x8": (124, 60),
    "hypar8x8": (128, 64),
    "test_surface_triang": (512, 0),
    "bunny170": (327, 0),
}


def sample(name):
    path = SAMPLES / f"{name}.obj"
    if not path.exists():
        pytest.skip(f"copy {name}.obj into examples/meshes/origamizer (see its README)")
    return og.load_mesh(path, triangulate_nonplanar=True)


@pytest.mark.parametrize("name", EXPECTED)
def test_sample_prepares_as_a_disk(name):
    prepared = sample(name)
    t = prepared.mesh.topology()
    assert t.is_disk  # every sample is an open surface with one boundary loop
    assert (t.faces, len(prepared.triangulated_faces)) == EXPECTED[name]
    assert math.isclose(sum(prepared.mesh.gaussian_curvature()), 2 * math.pi, abs_tol=1e-9)
    # A disk needs no cut, so the minimal cut adds no seams.
    assert og.cut_to_disk(prepared.mesh).cut_edges == []


@pytest.mark.parametrize("name", ["bell8x8", "hypar8x8", "test_surface_triang"])
def test_sample_maps_on_the_first_attempt(name):
    mapping = og.map_surface(sample(name).mesh)
    assert mapping.feasible, mapping.attempts[-1].message
    assert len(mapping.attempts) == 1  # the minimal cut, no remapping
    report = og.check_overlap(mapping.molecules, mapping.layout)
    assert report.valid and report.smallest_area > 0


@pytest.mark.parametrize("name", ["bell8x8", "hypar8x8"])
def test_grids_satisfy_t19(name):
    mapping = og.map_surface(sample(name).mesh, enforce_t19=True)
    assert mapping.feasible, mapping.attempts[-1].message
    conditions = og.crease_conditions(mapping.molecules, mapping.layout, mapping.theta, mapping.w)
    assert all(c.t19 for c in conditions)


def test_objective_b_on_the_test_surface_needs_remapping():
    # Phase 2 of objective B fails in HiGHS's QP solver with the plain angles
    # (a solve error, recorded as a failed attempt); the remapped angles work.
    mapping = og.map_surface(sample("test_surface_triang").mesh,
                             objective=og.WidthObjective.PAPER_PENALTY)
    assert mapping.feasible
    assert [a.remapped for a in mapping.attempts] == [False, True]


@pytest.mark.xfail(strict=True, reason="T1 and the convexity of vertex-tucking molecules (T8) "
                   "conflict at saddle vertices of bunny170 (docs/notes.tex, M2)")
def test_bunny170_maps():
    mapping = og.map_surface(sample("bunny170").mesh)
    assert mapping.feasible
