#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Matías Saavedra

"""Render the 2D layouts of the example meshes (milestone M2).

Run from the repository root:

    uv run --group plot python examples/plot_layout.py [--out build/plots] [--objective A|B]

For each mesh: the automatic cut, the rotation step, the translation step, and
the reconstruction of the molecule mesh (TVCG10 Sec. 3 and 6). Where no
nonnegative widths exist, the layout uses the widths of objective B's phase 1,
which solve T17 with some widths negative; where the rotation step itself is
infeasible, the panel says so. Faces are gray, edge-tucking
molecules blue (red where the width is negative, so the molecule is crossed),
vertex-tucking molecules yellow, and digons dark lines. Writes one PNG per
mesh and an overview, m2-layouts.png.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection, PolyCollection  # noqa: E402

import origamizer as og  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
MESHES = ROOT / "examples" / "meshes"
EXAMPLES = ["tetrahedron", "cube", "octahedron", "l_prism", "torus", "cylinder", "saddle"]
if (MESHES / "bunny.obj").exists():  # generated locally by fetch_bunny.py
    EXAMPLES.append("bunny")

FACE = (0.80, 0.81, 0.84)
ETM = (0.62, 0.78, 0.94)
CROSSED = (0.93, 0.35, 0.30)
VTM = (0.99, 0.90, 0.55)


def layout_of(name, objective):
    """Cut, molecule mesh, both mapping steps, and the reconstruction."""
    disk = og.cut_to_disk(og.load_mesh(MESHES / f"{name}.obj").mesh).disk
    mesh = og.build_molecule_mesh(disk)
    rotation = og.solve_rotation(mesh)
    if not rotation.feasible:
        return None, None, None, "rotation step infeasible:\n" + rotation.message.split(", largest")[0]
    translation = og.solve_translation(mesh, rotation.theta, objective=objective)
    status = "feasible"
    if not translation.feasible:
        translation = og.solve_translation(mesh, rotation.theta,
                                           objective=og.WidthObjective.PAPER_PENALTY)
        if not translation.w:
            raise SystemExit(f"{name}: translation step failed: {translation.message}")
        scale = og.mean_molecule_length(mesh)
        status = f"negative widths, total {translation.violation / scale:.2f} mean edges"
    layout = og.reconstruct_layout(disk, mesh, rotation.theta, translation.w)
    return mesh, translation.w, layout, status


def corner(layout, polygon, vertex):
    return layout.corners[polygon][layout.vertices[polygon].index(vertex)]


def draw(ax, name, objective):
    mesh, widths, layout, status = layout_of(name, objective)
    ax.set_title(f"{name}\n{status}", fontsize=9)
    ax.set_axis_off()
    if layout is None:
        print(f"{name}: {status}")
        return
    digons = range(mesh.faces, mesh.polygons)

    # Vertex-tucking molecules: the copies of each vertex, counterclockwise.
    vtms = [[corner(layout, slot.polygon, v) for slot in star.slots]
            for v, star in enumerate(mesh.stars)]
    # Edge-tucking molecules: the legs in the polygons before (P0) and after
    # (P1) the molecule, counterclockwise around its vertex i.
    etms, colors = [], []
    for star_vertex, star in enumerate(mesh.stars):
        for n, slot in enumerate(star.slots):
            molecule = mesh.molecules[slot.molecule]
            if molecule.kind == og.MoleculeKind.BOUNDARY or slot.sign != 1:
                continue
            p0 = slot.polygon
            p1 = star.slots[(n + 1) % len(star.slots)].polygon
            etms.append([corner(layout, p0, molecule.i), corner(layout, p0, molecule.j),
                         corner(layout, p1, molecule.j), corner(layout, p1, molecule.i)])
            colors.append(CROSSED if widths[slot.molecule] < 0 else ETM)
    faces = [layout.corners[p] for p in range(mesh.faces)]

    ax.add_collection(PolyCollection(vtms, facecolors=VTM, edgecolors=(0.6, 0.5, 0.2),
                                     linewidths=0.4))
    ax.add_collection(PolyCollection(etms, facecolors=colors, edgecolors=(0.3, 0.45, 0.6),
                                     linewidths=0.4, alpha=0.85))
    ax.add_collection(PolyCollection(faces, facecolors=FACE, edgecolors=(0.25, 0.25, 0.25),
                                     linewidths=0.7))
    ax.add_collection(LineCollection([layout.corners[p] for p in digons], colors="black",
                                     linewidths=1.6))
    xs = [x for polygon in layout.corners for x, _ in polygon]
    ys = [y for polygon in layout.corners for _, y in polygon]
    pad = 0.05 * max(max(xs) - min(xs), max(ys) - min(ys))
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(min(ys) - pad, max(ys) + pad)
    ax.set_aspect("equal")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=ROOT / "build" / "plots")
    parser.add_argument("--objective", choices=["A", "B"], default="A",
                        help="A: minimum total width (default); B: the paper's penalty")
    args = parser.parse_args()
    objective = {"A": og.WidthObjective.MINIMUM_WIDTH,
                 "B": og.WidthObjective.PAPER_PENALTY}[args.objective]
    args.out.mkdir(parents=True, exist_ok=True)

    for name in EXAMPLES:
        fig, ax = plt.subplots(figsize=(6, 6))
        draw(ax, name, objective)
        fig.savefig(args.out / f"layout_{name}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    columns = 4
    rows = -(-len(EXAMPLES) // columns)
    fig, axes = plt.subplots(rows, columns, figsize=(4 * columns, 4 * rows))
    for ax in axes.flat[len(EXAMPLES):]:
        ax.set_axis_off()
    for ax, name in zip(axes.flat, EXAMPLES):
        draw(ax, name, objective)
    fig.suptitle(f"M2 layouts, objective {args.objective}", fontsize=12)
    fig.savefig(args.out / "m2-layouts.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {len(EXAMPLES) + 1} images to {args.out}")


if __name__ == "__main__":
    main()
