# CLAUDE.md

Guidance for Claude Code in this repository. Read this file at the start of every session, check the Status log at the end, and update the log and the milestone checkboxes before the session ends.

## Project

- **Name:** Open Origamizer. Repository `MFSaavedra/open-origamizer` (https://github.com/MFSaavedra/open-origamizer), public from the first commit.
- **Goal:** a clean-room, open-source implementation of Tachi's Origamizer (folding one uncut sheet into a given polyhedral surface by tucking the excess material to one side), with exporters aimed at fabricating foldable surfaces. First fabrication target: rigid face tiles 3D printed onto fabric.
- **Why it has to be a reimplementation:** Tachi's Origamizer is proprietary and closed source, distributed as a Win32 binary (v0.47, 2008) licensed for non-commercial use with attribution. There is no public source to fork. The provably correct algorithm of Demaine and Tachi (SoCG 2017) was announced for a future version, but as of November 2026 the download page still offers only v0.47.
- **Prior art (checked 2026-11-23):** four GitHub repos mention "origamizer" and none implements tucking molecules. `KoyanoBunsho/mac_origamizer_open` (MIT) is a spanning-tree net unfolder and says it is not Origamizer; `cedix/origamizer` is an unrelated web app; `co-demos/flaplab` is an empty 2019 scaffold; `susiDX/Origamizer` is six unlicensed Python files from 2019 (code not read).

Tracks:
1. **Core:** reimplement the practical heuristic from Tachi, "Origamizing Polyhedral Surfaces" (IEEE TVCG 2010): tucking molecules, two-step mapping, edge splitting.
2. **Fabrication:** turn a layout into printable face tiles on fabric (spec below).
3. **Theory (stretch):** the SoCG 2017 provable construction and thickness-aware variants.

## Hard constraints

- Clean-room only. Implement from published papers. Never decompile, disassemble, or reverse-engineer Tachi's binary, and never commit it or files it produced. If the user runs it locally for comparison, that stays outside the repo and follows its license.
- Do not commit PDFs of IEEE or ASME papers (TVCG 2010, Yang and You 2024); link them instead. The SoCG 2017 paper is CC-BY 3.0 and is included in `docs/papers/` with attribution. `.gitignore` ignores `*.pdf` except for the files listed there.
- License: GPL-3.0-or-later (decided 2026-11-23, because CGAL, which is mostly GPL-3.0+, may be used freely). Start every source file with an `SPDX-License-Identifier: GPL-3.0-or-later` line and the line `Copyright (C) 2026 Matías Saavedra` (the year when the file is created). FOLD and Origami Simulator are MIT (checked 2026-11-23), so their code may be borrowed if their copyright notices are kept.
- Commits carry no Claude attribution. This is the user's standing instruction and overrides default attribution settings: author every commit as `Matías Saavedra <50891732+MFSaavedra@users.noreply.github.com>` with commit signing off (see Commands), and add no `Co-Authored-By` or `Claude-Session` lines.

## Working style

- Prose in docs, code comments, commit messages, and chat: precise and direct. Never use em-dashes.
- Long-form write-ups (derivations, design notes, anything paper-like) go in LaTeX under `docs/`. Markdown is fine for the README, this file, and short notes.
- Priority order: correct geometric formulation, then tests, then fabrication features. Keep applied tangents short.
- Every algorithmic function cites its source in the docstring (for example "TVCG10 Sec. 4" or "SoCG17 Sec. 5").
- Prefer explicit invariants and assertions over silent repairs. Every geometric routine gets a unit test.

## Decisions

Settled in session 1 (2026-11-23):
- **Repository:** `open-origamizer`, public from the first commit. The Python import name stays `origamizer`.
- **License:** GPL-3.0-or-later (see Hard constraints). Copyright holder: Matías Saavedra.
- **Language:** a C++ core with Python bindings (chosen over a pure-Python reference implementation).
- **C++ dependencies:** CGAL may be used freely; its GPL-3.0+ packages are compatible with the project license.
- **Contact with Tachi** (tachi at idea.c.u-tokyo.ac.jp): not contacted yet. Revisit before announcing the project anywhere, and record the outcome here. If he offers private source under restrictive terms, whoever reads it must not write this implementation.

Still open:
- **C++ stack details:** build system and packaging, bindings library, test frameworks, the split between C++ and Python, and the optimizer for M2 (decide that one after the TVCG10 notes).
- **Viewer location:** path or URL of the user's Origami Simulator fork (needed by M3).
- **Paper access in cloud sessions:** see Environment notes.

## Formats and tooling

- Input: OBJ, triangle or convex-polygon meshes. Keep polygon faces intact: trimesh triangulates them on load (tested with trimesh 5.2.2), so do not use it to read input meshes.
- Interchange: FOLD (https://github.com/edemaine/fold). Export a `creasePattern` frame (2D `vertices_coords`, `edges_vertices`, `edges_assignment`, `edges_foldAngle`, `faces_vertices`) and optionally a `foldedForm` frame with 3D coordinates. Checked against `doc/spec.md` on 2026-11-23: the `edges_assignment` letters are B (border), M (mountain), V (valley), F (flat), U (unassigned), C (cut or slit), and J (join); `edges_foldAngle` is in degrees in [-180, 180], positive for valley, negative for mountain, and zero for flat, unassigned, and border edges.
- SVG export compatible with Origami Simulator's importer (check its documentation for the color conventions).
- Fabrication: STL for the tiles plus a sidecar file with the fabric pause height.
- Viewer: the user's fork of Origami Simulator (upstream https://github.com/amandaghassaei/OrigamiSimulator), which imports FOLD and SVG. It does not detect self-intersection, so treat it as a visual sanity check, not a validator.

## Environment notes (Claude Code cloud sessions, as of 2026-11-23)

- The Default cloud environment reaches package registries and GitHub only. Paper hosts are blocked: origami.c.u-tokyo.ac.jp, tsg.ne.jp, www.tsg.ne.jp, doi.org, drops.dagstuhl.de, arxiv.org, scripts.mit.edu.
- The SoCG17 paper is in `docs/papers/`. arXiv papers (OriStitch, Ku and Demaine) can be read through the alphaXiv connector. The TVCG10 author PDF needs an upload or an environment that allows origami.c.u-tokyo.ac.jp.
- Sessions cannot create GitHub repositories or call account-wide GitHub APIs. Repository-scoped calls work once the repo is attached to the session.
- Installed: g++ 13.3, clang 18.1, CMake 3.28, Ninja 1.11, clang-format 18, Python 3.13, uv 0.11. Not installed: CGAL, Boost, GMP, MPFR, Eigen, any TeX engine.

## Planned layout

Written for the earlier pure-Python plan. Revise it for the C++ core once the stack details are settled.

```
.
├── CLAUDE.md
├── README.md
├── LICENSE
├── pyproject.toml
├── src/origamizer/
│   ├── mesh.py              # OBJ I/O, half-edge structure, checks, disk cutting
│   ├── tucking.py           # tucking-molecule parameters and validity conditions (TVCG10)
│   ├── layout.py            # two-step mapping, penalty optimization, edge splitting
│   ├── creases.py           # edge and vertex tucking molecules, mountain/valley assignment
│   ├── io_fold.py           # FOLD export and import
│   ├── io_svg.py            # SVG export for Origami Simulator
│   ├── metrics.py           # efficiency, minimum features, crease counts, layer estimates
│   ├── fabrication/tiles.py # beveled face tiles to STL, pause height
│   └── socg17/              # stretch: tuck proxy, squashing, streams, sites, Voronoi folding
├── tests/
├── examples/meshes/         # tetrahedron, cube, octahedron, saddle, low-poly bunny
└── docs/
    ├── tvcg2010-notes.md
    └── notes.tex
```

## Milestones

- [ ] **M0 Setup.** Tooling, license, README. Fetch the TVCG 2010 author PDF (link below) and write `docs/tvcg2010-notes.md` before any core code. Cover: the tucking-molecule parameters (an angle θ(i,j) and a distance w(i,j) per pair of adjacent polygons, with θ(j,i) = -θ(i,j), as in Tachi's slides; confirm the exact definitions in the paper), the equality and inequality conditions for a valid crease pattern, the two-step mapping (a rotational step, then a translational step; confirm the exact formulation), and edge splitting. Progress: license and README done (2026-11-23); tooling and notes pending.
- [ ] **M1 Mesh.** OBJ loading, manifold and orientability checks, convex facets (split nonconvex faces), cutting closed or higher-genus meshes into a topological disk (one cut path per handle), choice of tuck side (inside for closed surfaces). Test meshes: tetrahedron, cube, octahedron, a saddle patch, a low-poly bunny.
- [ ] **M2 Layout.** Two-step mapping with penalty-energy optimization in the translational step, non-overlap constraints, and edge splitting when infeasible. Optional term from OriStitch: a one-sided quadratic penalty enforcing a minimum separation between pieces, weighted by λ.
- [ ] **M3 Crease pattern.** Edge-tucking and vertex-tucking molecules, mountain and valley assignment, FOLD and SVG export. Acceptance: the tetrahedron and the cube load in the viewer and fold toward the target with tucks inside.
- [ ] **M4 Metrics.** Area efficiency (surface area over paper area), crease count, smallest crease segment and smallest gap, estimated maximum layers per tuck. Reference point: the 374-triangle Stanford bunny folds from paper at 22.3% area efficiency, about 2:1 linear scale (SoCG17 Fig. 1).
- [ ] **M5 Fabrication.** Tile generator (spec below), minimum-feature and bed-size checks, a printed test on a cube and a saddle.
- [ ] **M6 Stretch.** SoCG17 construction on small meshes; thickness-aware objectives (see Research questions); notes in LaTeX.

## Algorithm reference: SoCG 2017 (Demaine and Tachi)

Definitions
- **Q:** an embedded polyhedral manifold with strictly convex facets, homeomorphic to a disk, with a clean side and a tuck side. It must not touch itself on the tuck side. Preprocessing: double the faces of nonorientable or nonmanifold inputs, cut handles, split nonconvex facets.
- **ε-extra folding** of a convex paper P into Q: the image contains Q and otherwise lies in the tuck-side ε offset of Q (Hausdorff distance at most ε).
- **Seamless:** each facet's clean side is covered by one crease-pattern face on the outermost layer. The construction gives more: every facet is a single uncreased face with no additional layers, even on the tuck side.
- **Watertight:** a closed curve C on P (here the boundary of P) folds to within Fréchet distance ε of the boundary of Q, and the folded interior of C covers Q.
- **Gauss-Bonnet obstruction:** a watertight folding cannot realize a negative-curvature vertex exactly, hence the small extra structure on the tuck side.

Pipeline
1. **Sec. 3 Tuck proxy T (a polyhedral waffle).** Inset the edges of Q on the tuck side by ε' (like the start of a 3D straight skeleton). The inset bisects every dihedral angle, splitting reflex ones into two convex ones, so the walls lie in the bisector planes. Around each vertex, connect the offsets by a cycle on a small sphere, then overlay a regular tetrahedron and triangulate so every pocket face has spherical perimeter at most 328.5 degrees (intrinsically convex pockets). Wall height is at most ε, with 0 < ε' < ε.
2. **Sec. 4 Local squashing.** For each waffle pocket, vertex-unfold the walls into the plane and place rays in the angular gaps so consecutive rays sandwich each wall and span less than 180 degrees. Intersect with a circle (vertex pockets) or with the facet offset whose corners are smoothed by quadratic Bezier curves (facet pockets). The result is a strictly convex cell; bottom angles only increase; cell perpendiculars run counterclockwise.
3. **Sec. 5 Stream placement.** Tutte-embed the waffle dual inside a convex polygon P' with the boundary nodes on its boundary. Wall quadrangles become rivers (thickened by the shared edge length); wall triangles become brooks (zero width). Around each node, nested disks R_x, S_x, T_x separate, bend, and twist the streams (k tracks for degree k). Fillet corners so streams are C1 with clearance δ > 0, and scale the embedding up until the disks are local.
4. **Sec. 6 Site placement.** Point, segment, and polygon sites; a generalized Voronoi diagram in the geodesic vertex metric (shortest-path distance to the nearest site vertex, with sites as obstacles). Requirements: (1) the Voronoi diagram contracts to the waffle graph of T; (2) every Voronoi edge is a straight segment that mirror-reflects its site vertices, giving paired subcells; (3) subcell angles are at least the corresponding wall bottom angles; (4) every point of P is within ε' of a site vertex. Protective squashes: pumpkins at brook nodes, jack-o'-lanterns at placed facets, calabashes along brooks, cucumbers along rivers. Order: node sites, then stream sites (bisect gaps as needed), then bank sites until every point is within ε'. P is the convex hull of the paired subcells.
5. **Sec. 7 Voronoi folding.** Fold P into an abstract waffle, which need not embed in 3D: paired subcells double-cover walls, unpaired subcells single-cover boundary walls. Contract non-realizing edges by folding wall triangles in half; unify parallel edges by crimping to the smallest angle; merge each realizing path by crimping bottom angles to those of T and stacking walls. The result embeds in T.

Caveats: the construction guarantees feasibility, not efficiency. The paper does not optimize within its family of foldings and expects that to be hard, so a direct implementation would likely give large sheets with fine creases. The original heuristic sometimes fails to find a feasible layout; that is the gap this construction closes.

## Fabrication track: rigid tiles on fabric

Rationale
- Facets never carry stacked layers, so they can be thick rigid tiles. All the difficulty is in the tucks, so the tucks are bare fabric.
- Out of scope: printing the whole crease pattern as one rigid thick sheet. The tucks stack zero-thickness layers inside an ε band, Origamizer certifies a final folded state but no folding motion, and rigid thick panels lose mobility (Yang and You: a uniform-thickness Miura unit with revolute joints has m = 6(9 - 1) - 12*5 = -12 and needs slits to move).

Tile spec, version 1
- One tile per facet at its laid-out 2D position, thickness T (parameter).
- Print with the clean side on the bed; pause after the first layers, lay tulle or lycra flat, resume. Output the pause height.
- Bevel tile edges where neighbours would collide when folded. For a tuck-side dihedral angle φ < π, cut the edge face at (π - φ)/2 from the tile normal so it lies along the dihedral bisector (45 degrees for a cube). That bisector is where the SoCG17 tuck proxy puts the walls, so the fabric fin lands between the bevels. Offset each bevel by one fabric thickness (total gap about two fabric thicknesses). Measure bevels from the fold line in the fabric plane; the thin skin under the fabric opens slightly at convex edges.
- Tuck regions stay bare fabric. Optional mode, following OriStitch: mark vertex-tuck regions for cutting out, which lowers closing forces at the cost of watertightness.
- Checks: warn when a facet edge is below the minimum (start at 8 to 9 mm; OriStitch's embroidery floor is 8.4 mm) or when the layout exceeds the bed (default 250 x 250 mm).
- Closing: by hand or with drawstrings at first. Boundary faces stay floppy; plan a printed rim (OriStitch sewed one on).

Scale: on a 250 mm bed at 22.3% efficiency, a 374-facet model gets about 0.37 cm2 per facet (about 9 mm across) and edge tucks only a few millimetres deep. Keep first prints to tens of faces.

Lessons from OriStitch (Origamizer-based, fabric with heat-shrink embroidery)
- Modified only the translational-mapping phase, adding a minimum-separation penalty.
- Four hinge types chosen by gap width; the narrowest gaps are plain creases.
- Laser-cut the material between vertices, so no vertex-tucking molecules were needed. The uncut layout also folds but needs more force.
- Saddle-like regions (low mean curvature, nonzero Gaussian curvature) push faces far apart, giving cycles of wide hinges with heavy crimping.
- Converted 26 of 28 benchmark meshes; both failures came from remeshing to the minimum resolution.
- Simultaneous actuation (oven) made hinges pull on each other; sequential heat-gun actuation worked.

## Research questions (theory track)

- **Layer count.** Minimize the maximum number of layers in any wall subject to watertightness. In SoCG17 Sec. 7 each Voronoi edge along a realizing path contributes a double-covered wall and contracted non-realizing edges add more, so layer count is governed by site placement.
- **Clearance.** Certify vertex-pocket clearance for a material thickness t.
- **Reachability.** Is the target state reachable with rigid faces and flexible tucks? With all creases actuated at once? OriStitch notes that valid folding motions are often not known to exist.
- **Related prior art.** Ku and Demaine's offset-crease thickening (convex faces, no crease touching another crease's interior, bounds on thickness; connectivity of the flat and folded states only shown numerically at IASS 2016). Yang and You's slit method for uniform-thickness Miura-ori (one DoF from merged 8R loops; compact folding needs 2b cos α ≥ c).

## References

Core
- Demaine, Tachi. Origamizer: A Practical Algorithm for Folding Any Polyhedron. SoCG 2017. https://doi.org/10.4230/LIPIcs.SoCG.2017.34 (CC-BY 3.0; copy in `docs/papers/`)
- Tachi. Origamizing Polyhedral Surfaces. IEEE TVCG 16(2):298-311, 2010. Author PDF: https://origami.c.u-tokyo.ac.jp/~tachi/cg/OrigamizingPolyhedralSurfaces_tachi_IEEETVCG.pdf
- Tachi. 3D Origami Design based on Tucking Molecule. Origami4, pp. 259-272, 2009. http://www.tsg.ne.jp/TT/cg/3DOrigamiDesign_tachi_4OSME.pdf
- Tachi's software page and license: https://tsg.ne.jp/TT/software/

Fabrication
- Chang et al. OriStitch: A Machine Embroidery Workflow to Turn Existing Fabrics into Self-Folding 3D Textiles. arXiv:2412.02891. https://arxiv.org/abs/2412.02891
- Narumi et al. Inkjet 4D Print: Self-folding Tessellated Origami Objects by Inkjet UV Printing. ACM TOG 42(4), 2023. https://doi.org/10.1145/3592409
- Suto et al. Crane (Grasshopper platform; thickening with hinges). ACM TOCHI 2023. https://pp.ics.keio.ac.jp/?p=171
- Segall, Ren, Vaxman, Sorkine-Hornung. Fabric Tessellation: Realizing Freeform Surfaces by Smocking. ACM TOG 43(4), 2024. https://doi.org/10.1145/3658151
- BYU CMR printable thick-origami examples: https://www.compliantmechanisms.byu.edu/maker-resources

Thickness
- Yang, You. Folding Miura-Ori of Uniform Thickness With One Degree-of-Freedom. ASME JMR 16(11):111004, 2024. https://doi.org/10.1115/1.4065004
- Ku, Demaine. Folding Flat Crease Patterns With Thick Materials. ASME JMR 8(3), 2016. https://arxiv.org/abs/1601.05747
- Ku, Demaine. Rigid folding analysis of offset crease thick folding. IASS 2016. https://scripts.mit.edu/~jasonku/pdf/RIGID_IASS2016.pdf
- Tuck-Folding Method for the Design of Origami Structures with Creased Folds (Springer chapter; small supplementary ZIP, license unclear). https://link.springer.com/chapter/10.1007/978-3-319-91866-2_4

Tools
- FOLD format: https://github.com/edemaine/fold
- Origami Simulator (upstream): https://github.com/amandaghassaei/OrigamiSimulator

## Commands

Not set up yet. Once tooling exists, list the exact commands here: install, test, run an example, and export FOLD, SVG, and STL.

Git identity. A fresh cloud container starts with a Claude identity and commit signing on, so run these in the repo after every clone:

```bash
git config user.name "Matías Saavedra"
git config user.email "50891732+MFSaavedra@users.noreply.github.com"
git config commit.gpgsign false
```

## Status log

- 2026-11-23: Project scoped in a claude.ai conversation. No code yet; repo not created. Next: settle the open decisions, then M0.
- 2026-11-23 (session 1, Claude Code): Created the public repo with this file, README, GPL-3.0 license, `.gitignore` guardrails, and the SoCG17 paper. Settled: name, visibility, license, language (C++ core with Python bindings), CGAL allowed, commit authorship. Tachi not contacted yet. Next: settle the C++ stack details, then M0.
