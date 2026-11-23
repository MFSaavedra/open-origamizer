#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Matías Saavedra

"""Render the example meshes with their automatic cuts.

Run from the repository root:

    uv run --group plot python examples/plot_cut.py [--out build/plots]

Faces are light gray, the original boundary is blue, and cut edges are red.
Writes one PNG per mesh and an overview, m1-cuts.png. matplotlib does not sort
depth exactly, so edges behind a face can show through; this is a visual
check, not a renderer.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection  # noqa: E402

import origamizer as og  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
MESHES = ROOT / "examples" / "meshes"
EXAMPLES = ["tetrahedron", "cube", "octahedron", "l_prism", "torus", "cylinder", "saddle"]
if (MESHES / "bunny.obj").exists():  # generated locally by fetch_bunny.py
    EXAMPLES.append("bunny")
# Elevation, azimuth, and vertical axis (the Stanford bunny is stored with y up).
VIEWS = {"saddle": (35, -60, "z"), "torus": (45, -50, "z"), "cylinder": (25, -60, "z"),
         "bunny": (15, 25, "y")}


def draw(ax, name):
    mesh = og.load_mesh(MESHES / f"{name}.obj").mesh
    cut = og.cut_to_disk(mesh)
    points = mesh.points

    polygons = [[points[i] for i in face] for face in mesh.faces]
    ax.add_collection3d(Poly3DCollection(polygons, facecolor=(0.82, 0.84, 0.9, 0.5),
                                         edgecolor=(0.35, 0.35, 0.35, 0.7), linewidths=0.6))
    boundary = [(points[a], points[b]) for loop in mesh.boundary_loops()
                for a, b in zip(loop, loop[1:] + loop[:1])]
    if boundary:
        ax.add_collection3d(Line3DCollection(boundary, colors="tab:blue", linewidths=2.2))
    cuts = [(points[a], points[b]) for a, b in cut.cut_edges]
    if cuts:
        ax.add_collection3d(Line3DCollection(cuts, colors="tab:red", linewidths=3.0))

    lows = [min(p[k] for p in points) for k in range(3)]
    highs = [max(p[k] for p in points) for k in range(3)]
    half = max(h - lo for h, lo in zip(highs, lows)) / 2
    centers = [(h + lo) / 2 for h, lo in zip(highs, lows)]
    ax.set_xlim(centers[0] - half, centers[0] + half)
    ax.set_ylim(centers[1] - half, centers[1] + half)
    ax.set_zlim(centers[2] - half, centers[2] + half)
    ax.set_box_aspect((1, 1, 1))
    elevation, azimuth, vertical = VIEWS.get(name, (22, -55, "z"))
    ax.view_init(elevation, azimuth, vertical_axis=vertical)
    ax.set_axis_off()
    t = mesh.topology()
    ax.set_title(f"{name}\nchi {t.euler_characteristic}, genus {t.genus}, "
                 f"{len(cut.cut_edges)} cut edges", fontsize=9)


def main():
    parser = argparse.ArgumentParser(description="Render the example meshes with their cuts.")
    parser.add_argument("--out", type=Path, default=ROOT / "build" / "plots")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    for name in EXAMPLES:
        fig = plt.figure(figsize=(4, 4))
        draw(fig.add_subplot(projection="3d"), name)
        fig.savefig(args.out / f"{name}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    fig = plt.figure(figsize=(16, 8.5))
    for k, name in enumerate(EXAMPLES):
        draw(fig.add_subplot(2, 4, k + 1, projection="3d"), name)
    fig.suptitle("Automatic cuts (red) and original boundaries (blue)", fontsize=13)
    fig.savefig(args.out / "m1-cuts.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {len(EXAMPLES) + 1} images to {args.out}")


if __name__ == "__main__":
    main()
