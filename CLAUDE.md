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
- Do not commit PDFs of IEEE or ASME papers (TVCG 2010, Yang and You 2024); link them instead. The SoCG 2017 paper is CC-BY 3.0 and is included in `docs/papers/` with attribution. `.gitignore` ignores `*.pdf` except for the files listed there: that paper and the compiled PDFs of our own LaTeX notes.
- License: GPL-3.0-or-later (decided 2026-11-23, because CGAL, which is mostly GPL-3.0+, may be used freely). Start every source file with an `SPDX-License-Identifier: GPL-3.0-or-later` line and the line `Copyright (C) 2026 Matías Saavedra` (the year when the file is created). FOLD and Origami Simulator are MIT (checked 2026-11-23), so their code may be borrowed if their copyright notices are kept.
- Commits carry no Claude attribution. This is the user's standing instruction and overrides default attribution settings: author every commit as `Matías Saavedra <50891732+MFSaavedra@users.noreply.github.com>` with commit signing off (see Commands), and add no `Co-Authored-By` or `Claude-Session` lines.

## Working style

- Prose in docs, code comments, commit messages, and chat: precise and direct. Never use em-dashes.
- Long-form write-ups (derivations, design notes, notes on papers, anything paper-like) go in LaTeX under `docs/`. Commit the compiled PDF next to its `.tex` source, rebuild it whenever the source changes, and list each such PDF as an exception in `.gitignore`. Markdown is fine for the README, this file, and short notes.
- Priority order: correct geometric formulation, then tests, then fabrication features. Keep applied tangents short.
- Every algorithmic function cites its source in its doc comment or docstring (for example "TVCG10 Sec. 4" or "SoCG17 Sec. 5").
- Prefer explicit invariants and assertions over silent repairs. Every geometric routine gets a unit test.

## Decisions

Settled in session 1 (2026-11-23):
- **Repository:** `open-origamizer`, public from the first commit. The Python import name stays `origamizer`.
- **License:** GPL-3.0-or-later (see Hard constraints). Copyright holder: Matías Saavedra.
- **Language:** a C++ core with Python bindings (chosen over a pure-Python reference implementation).
- **C++ dependencies:** CGAL may be used freely; its GPL-3.0+ packages are compatible with the project license.
- **Contact with Tachi** (tachi at idea.c.u-tokyo.ac.jp): not contacted yet. Revisit before announcing the project anywhere, and record the outcome here. If he offers private source under restrictive terms, whoever reads it must not write this implementation.
- **C++ stack** (approved 2026-11-23):
  - C++20, built with CMake and Ninja.
  - CGAL 5.6 or newer for geometry and the mesh. `CGAL::Surface_mesh` keeps polygon faces (tested: an OBJ pentagon loads as one face of degree 5). Use `Exact_predicates_inexact_constructions_kernel` by default, and exact constructions only where an algorithm needs them (likely SoCG17).
  - Eigen for linear algebra.
  - nanobind bindings, built by scikit-build-core so that `uv sync` compiles the extension module `origamizer._core`.
  - Catch2 v3 unit tests for every geometric routine, run under ASan and UBSan in CI; pytest for the bindings and end-to-end runs (OBJ in, FOLD out).
  - Split: algorithms in C++ (mesh, tucking, layout, creases, metrics); Python for the FOLD, SVG, and STL writers, fabrication, and the CLI.
  - Invariants: always-on checks that throw `origamizer::InvariantError`, which nanobind turns into a Python exception. Never plain `assert`, which release builds (wheels included) compile out.

- **Git workflow** (2026-11-23): direct pushes to `main` for now.
- **Notes format** (2026-11-23): LaTeX for all long-form notes, notes on papers included, with the compiled PDF committed next to the source. Design decisions beyond the papers go in `docs/notes.tex`, one section per milestone.
- **Clean side** (2026-11-23): face normals point to the clean side and tucks go behind. Closed meshes are oriented outward, so tucks go inside; meshes with boundary keep their orientation (TVCG10 Fig. 20).
- **Test data** (2026-11-23): the Stanford bunny is generated locally by `examples/meshes/fetch_bunny.py` and never committed, because the Stanford terms forbid commercial use without permission.

- **M2 solver** (2026-11-24): HiGHS v1.15.1 (MIT), fetched and built from source by CMake at that pinned tag: static, not sanitized, its warnings not fatal. The user first preferred CGAL's exact LP/QP solver, but on synthetic problems shaped like the translation step it took 0.05 s at 30 variables and over 10 minutes at 300 (exact arithmetic; with floating point it aborts), while HiGHS took 0.01 to 0.02 s up to 600 variables. HiGHS also beats the paper's projected gradient descent: it proves infeasibility and solves both objectives directly. v1.15.1 rather than v1.16.0 because v1.16.0 was 13 days old.
- **M2 objectives** (2026-11-24): implement both, selectable: (A) minimize the total tuck width subject to a minimum width, an LP that turns OriStitch's minimum separation into a hard bound; (B) the paper's penalty (TVCG10 Sec. 6.3) as quadratic programs. Exact definitions and default parameters are in the M2 section of `docs/notes.tex`.
- **M2 defaults** (2026-11-24, accepted by the user as proposed): boundary molecules follow our Gauss-Bonnet-consistent reading (TVCG10 notes, Sec. 4); convex faces stay whole, with triangulation as a fallback if the cube's quadrilaterals cannot be laid out; 2D conditions first (T1 to T9, T19, edge splitting), the tuck proxy last.
- **Boundary digons on by default** (2026-11-24; changes the proposed default, which had them optional): T1 with T8 needs at least three molecules around every vertex, and every example has boundary vertices with one or two faces (all 40 on the torus, 44 of 114 on the bunny), so without the TVCG10 Sec. 5.1 digons the rotation step is infeasible on every example. Derivation and counts in `docs/notes.tex`, M2. `MoleculeOptions::boundary_digons` can turn them off.
- **Objective B, phase 1 is a linear program** (2026-11-24; changes the accepted definition): it minimizes the total negative width instead of the sum of squared negative widths. HiGHS's QP solver failed on the squared version, whose widths are free variables with zero curvature (a wrong optimum on the cube, solve errors on three examples, an internal assertion on the tetrahedron). Which layouts count as feasible is unchanged; only the widths reported for an infeasible layout differ (least total violation instead of least squared violation).
- **Solver accuracy** (2026-11-24): every HiGHS optimum is polished (the KKT system of its active set solved exactly with a sparse QR) and must pass an independent KKT check, or the solve fails with a message. HiGHS's QP solver regularizes the Hessian and stops near 1e-7, and that noise broke exact structure: with it, objective A called the cube feasible with widths near 1e8. Details in `docs/notes.tex`, M2.

Still open:
- **Viewer location:** path or URL of the user's Origami Simulator fork (needed by M3).
- **Open questions in the TVCG10 notes:** seven, listed in their last section. Settle each before the milestone that needs it.
- **Cut quality:** the automatic cut is only guaranteed valid. Compare cut strategies once M2 can measure layouts (see the M1 remark in `docs/notes.tex`).

## Formats and tooling

- Input: OBJ, triangle or convex-polygon meshes. Keep polygon faces intact: trimesh triangulates them on load (tested with trimesh 5.2.2), so do not use it to read input meshes.
- Interchange: FOLD (https://github.com/edemaine/fold). Export a `creasePattern` frame (2D `vertices_coords`, `edges_vertices`, `edges_assignment`, `edges_foldAngle`, `faces_vertices`) and optionally a `foldedForm` frame with 3D coordinates. Checked against `doc/spec.md` on 2026-11-23: the `edges_assignment` letters are B (border), M (mountain), V (valley), F (flat), U (unassigned), C (cut or slit), and J (join); `edges_foldAngle` is in degrees in [-180, 180], positive for valley, negative for mountain, and zero for flat, unassigned, and border edges.
- SVG export compatible with Origami Simulator's importer (check its documentation for the color conventions).
- Fabrication: STL for the tiles plus a sidecar file with the fabric pause height.
- Viewer: the user's fork of Origami Simulator (upstream https://github.com/amandaghassaei/OrigamiSimulator), which imports FOLD and SVG. It does not detect self-intersection, so treat it as a visual sanity check, not a validator. Plan for M3 (agreed 2026-11-23): export FOLD with each crease's target fold angle taken from our folded state; read the simulator's FOLD importer source (MIT) first to confirm the fields it reads and how it handles polygon faces; then drive it headlessly with Playwright (Chromium is preinstalled in cloud sessions) to screenshot fold stages. Before M3, render the mesh, cut, and layout with matplotlib.

## Environment notes (Claude Code cloud sessions, as of 2026-11-23)

- Sessions run in the user's Default cloud environment. It now has unrestricted network access (it was limited to package registries and GitHub when session 1 started), so the TVCG10 author PDF, www.tsg.ne.jp, and arXiv all download. Use `https://` URLs: plain HTTP fails because the container reaches the internet only through the HTTPS proxy.
- The SoCG17 paper is in `docs/papers/`. The alphaXiv connector also serves arXiv papers.
- CMake fetches HiGHS from GitHub at configure time (about 20 s) and compiles it with the project (about 100 s more on a clean build); `uv sync` does the same in its own build directory. On GitHub Actions, compiling HiGHS adds 1.5 to 3 minutes to each job on a cold cache; since 2026-11-24 the workflow compiles through ccache, kept between runs with `actions/cache` (first run 45307715973, cold). The apt step varies far more, from 16 s to 6 minutes between runs (mirror speed), and often dominates the job time. An earlier note here blamed HiGHS alone for jobs of about 8 minutes; that was wrong.
- Every session starts from a fresh container, and neither the C++ dependencies nor TeX are preinstalled. Install them with the apt commands in Commands (about 10 s and 50 s). A setup script could do this automatically, but it would need a separate cloud environment for this repo; none is set up.
- GitHub refuses pushes from these sessions that add or change files under `.github/workflows/`, because the session token lacks the `workflow` scope. The user adds such files on github.com.
- Sessions cannot create GitHub repositories or call account-wide GitHub APIs. Repository-scoped calls work once the repo is attached to the session.
- Preinstalled: g++ 13.3, clang 18.1, CMake 3.28, Ninja 1.11, clang-format 18, Python 3.13, uv 0.11. Not preinstalled: CGAL, Boost, GMP, MPFR, Eigen, Catch2, any TeX engine.

## Planned layout

M0 and M1 created the files marked as existing below; the rest is planned.

```
.
├── CLAUDE.md
├── README.md
├── LICENSE
├── .github/workflows/ci.yml # CI: C++ tests under sanitizers (g++, clang++), Python tests
├── CMakeLists.txt
├── pyproject.toml           # scikit-build-core + nanobind, managed with uv
├── include/origamizer/      # public C++ headers
│   ├── check.hpp            # (exists) always-on invariant checks (throw InvariantError)
│   ├── mesh.hpp             # (exists) OBJ input, checks, convex splitting, topology, clean side
│   ├── cut.hpp              # (exists) cutting into a topological disk
│   ├── molecules.hpp        # (exists) molecule mesh: edge, digon, and boundary molecules, vertex stars
│   ├── layout.hpp           # (exists: rotation, translation, reconstruction) remapping, edge splitting
│   ├── tucking.hpp          # tuck proxy and 3D conditions (TVCG10 Sec. 4.3)
│   ├── creases.hpp          # edge and vertex tucking molecules, mountain/valley assignment
│   ├── metrics.hpp          # efficiency, minimum features, crease counts, layer estimates
│   └── socg17/              # stretch: tuck proxy, squashing, streams, sites, Voronoi folding
├── src/                     # C++ sources matching the headers
│   ├── detail/              # (exists) private helpers: union-find, the HiGHS wrapper
│   └── bindings/            # (exists) nanobind module origamizer._core
├── python/origamizer/       # Python package
│   ├── io_fold.py           # FOLD export and import
│   ├── io_svg.py            # SVG export for Origami Simulator
│   ├── fabrication/tiles.py # beveled face tiles to STL, pause height
│   └── cli.py
├── tests/
│   ├── cpp/                 # Catch2 unit tests
│   └── python/              # pytest: bindings and end-to-end runs
├── examples/
│   ├── meshes/              # (exists) generate.py and its OBJ files; fetch_bunny.py (bunny.obj ignored)
│   └── plot_cut.py          # (exists) renders each example with its cut
└── docs/
    ├── papers/              # openly licensed papers only
    ├── tvcg2010-notes.tex   # (exists) and its compiled .pdf
    └── notes.tex            # (exists) design notes, one section per milestone, and its .pdf
```

## Milestones

- [x] **M0 Setup.** Tooling, license, README. Fetch the TVCG 2010 author PDF (link below) and write `docs/tvcg2010-notes.tex` (with its compiled PDF) before any core code. Cover: the tucking-molecule parameters (an angle θ(i,j) and a distance w(i,j) per pair of adjacent polygons, with θ(j,i) = -θ(i,j), as in Tachi's slides; confirm the exact definitions in the paper), the equality and inequality conditions for a valid crease pattern, the two-step mapping (a rotational step, then a translational step; confirm the exact formulation), and edge splitting. Done 2026-11-23: build skeleton, tests, and `docs/tvcg2010-notes.tex` with its PDF. CI runs on GitHub Actions; its first run passed the same day (C++ tests under ASan and UBSan with g++ and clang++, Python tests).
- [x] **M1 Mesh.** OBJ loading, manifold and orientability checks, convex facets (split nonconvex faces), cutting closed or higher-genus meshes into a topological disk (a cut graph with 2g loops for genus g, plus paths joining the boundary loops; "one cut path per handle" was wrong), choice of tuck side (inside for closed surfaces). Test meshes: tetrahedron, cube, octahedron, a saddle patch, a low-poly bunny. Done 2026-11-23: `mesh.hpp` and `cut.hpp` with C++ and Python tests, example meshes plus a torus, an open cylinder, and an L prism, the bunny via `fetch_bunny.py`, renders via `plot_cut.py`, and design notes in `docs/notes.tex`. Cut quality is deferred to M2.
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
- Tachi. 3D Origami Design based on Tucking Molecule. Origami4, pp. 259-272, 2009. https://www.tsg.ne.jp/TT/cg/3DOrigamiDesign_tachi_4OSME.pdf
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

Verified in a cloud session on 2026-11-23. Commands to run an example and to export FOLD, SVG, and STL come with M1 to M5.

System packages (Ubuntu 24.04 cloud container, run as root). C++ dependencies, about 10 s: CGAL 5.6 with Boost, GMP, and MPFR, plus Eigen 3.4 and Catch2 3.4.

```bash
apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends libcgal-dev libeigen3-dev catch2
```

LaTeX for the notes, about 50 s:

```bash
DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends texlive-latex-recommended texlive-latex-extra texlive-fonts-recommended texlive-pictures texlive-science lmodern latexmk
```

C++ tests in Debug with AddressSanitizer, UBSan, and warnings as errors (the CI configuration):

```bash
cmake -S . -B build/cpp-asan -G Ninja -DCMAKE_BUILD_TYPE=Debug -DORIGAMIZER_SANITIZE=ON -DCMAKE_COMPILE_WARNING_AS_ERROR=ON
cmake --build build/cpp-asan
ctest --test-dir build/cpp-asan --output-on-failure
```

Python package and its tests. `uv sync` compiles `origamizer._core` and rebuilds it when C++ sources or build files change:

```bash
uv sync
uv run pytest
```

Notes. Rebuild the PDF after every change to a `.tex` file and commit both. The build is reproducible, so an unchanged source gives an identical PDF:

```bash
latexmk -pdf -cd docs/tvcg2010-notes.tex
latexmk -pdf -cd docs/notes.tex
```

Example meshes: regenerate the committed ones, fetch and simplify the bunny (needs network; output not committed), and render every mesh with its cut into `build/plots/`:

```bash
python3 examples/meshes/generate.py
uv run python examples/meshes/fetch_bunny.py
uv run --group plot python examples/plot_cut.py
```

CI status from a session (repository-scoped GitHub calls work; check that the test counts in the log are nonzero, not just that the run is green):

```bash
gh run list -R MFSaavedra/open-origamizer --limit 5
gh run view <run-id> -R MFSaavedra/open-origamizer --log
```

Git identity. A fresh cloud container starts with a Claude identity and commit signing on, so run these in the repo after every clone:

```bash
git config user.name "Matías Saavedra"
git config user.email "50891732+MFSaavedra@users.noreply.github.com"
git config commit.gpgsign false
```

## Status log

- 2026-11-23: Project scoped in a claude.ai conversation. No code yet; repo not created. Next: settle the open decisions, then M0.
- 2026-11-23 (session 1, Claude Code): Created the public repo with this file, README, GPL-3.0 license, `.gitignore` guardrails, and the SoCG17 paper. Settled: name, visibility, license, language (C++ core with Python bindings), CGAL allowed, commit authorship. Tachi not contacted yet. Later the same session: the user approved the C++ stack (see Decisions). The Default environment now has unrestricted network access, so the TVCG10 PDF is reachable. The C++ dependencies install with apt in about 10 s, and a CGAL 5.6 test program compiled and kept an OBJ pentagon as one face. Next: M0 (tooling, then the TVCG10 notes).
- 2026-11-23 (session 1, continued): Decided direct pushes to `main` and LaTeX notes with committed PDFs. M0 done: build skeleton with always-on invariant checks (7 Catch2 tests passing under ASan and UBSan with g++ and clang++, 2 pytest tests), and the TVCG10 notes (11 pages, PDF committed; seven open questions). The CI workflow could not be pushed to `.github/workflows/` (no `workflow` scope), so it is parked at `ci/ci.yml` for the user to move on github.com. Next: M1 (mesh).
- 2026-11-23 (session 1, continued): The user moved the workflow to `.github/workflows/ci.yml`. First CI run (45225589519) passed all three jobs: 7/7 Catch2 tests under ASan and UBSan with g++ and with clang++, 2/2 pytest. Only warnings: GitHub's Node.js 20 deprecation notice for `actions/checkout@v4` and `astral-sh/setup-uv@v6`; v7 of both runs on Node 24 and needs no new inputs, and the user decides whether to edit the workflow on github.com. Next: M1 (mesh).
- 2026-11-23 (session 1, continued): M1 done. Mesh input, checks, convex splitting, topology, clean side, and disk cutting in C++ with Python bindings: 48 Catch2 tests (ASan and UBSan) and 24 pytest tests, the bunny test only when `bunny.obj` exists. The invariants caught two bugs during development, both fixed: one-edge seed slits were glued back, and cut copies of a vertex tripped the coincident-vertex check. The Stanford archive's coarse bunnies are not manifold, so the full one is simplified to 372 triangles; it has 5 scan holes and needs 45 cut edges. Origami Simulator plan recorded under Formats and tooling. Next: M2 (layout), starting with the M2 optimizer decision and TVCG10 notes Open questions 1 and 5.
- 2026-11-24 (session 2): M2 started. Benchmarked CGAL's exact LP/QP solver against HiGHS on synthetic layout-sized problems (CGAL over 10 minutes at 300 variables, HiGHS 0.01 s); the user chose HiGHS v1.15.1 (MIT), now fetched and built by CMake behind a private wrapper. Decided to implement both translation objectives (minimum width; the paper's penalty) and accepted the proposed defaults; the boundary-digon default was then changed with evidence (T1 with T8 needs at least three molecules per vertex; see Decisions). Done: molecule mesh (`molecules.hpp`) and the rotation step (`layout.hpp`), feasible on all seven committed examples with the Gauss-Bonnet boundary turning verified on every solution: 61 Catch2 tests. Next: translation step (objectives A and B, and the rank of C_w), polygon placement with the overlap check and renders, then T19, edge splitting, remapping, and the tuck proxy.
- 2026-11-24 (session 2, continued): Translation step done, with both objectives, and `reconstruct_layout`, which places every polygon and checks every side of every vertex-tucking molecule against T2 and T3; on all seven examples the layout closes to 2e-13 of the mean edge length, which verifies T17 independently. Found and fixed: HiGHS treats matrix entries up to 1e-9 as zero with a warning (the wrapper drops them), its QP optima carry 1e-8 noise that broke exact zeros and dependencies (the wrapper now polishes every optimum and checks the KKT conditions), and the squared phase 1 of objective B made its QP solver fail (now a linear program; see Decisions). C_w has rank 2V - 1 on four examples, consistent each time (TVCG10 notes, Open question 5, now answered). Result: with the rotation step's angles, nonnegative widths exist on five examples but not on the cube (total negative width 11.46 mean edge lengths) or the L prism (8.16), the TVCG10 Sec. 6.4 failure. 71 Catch2 tests under ASan and UBSan with g++, and with clang in Release. Next: Python bindings and renders of the layouts, then remapping (TVCG10 Sec. 6.5), which needs the rotation step to minimize the distance to reference angles; then T19, edge splitting, and the tuck proxy.
