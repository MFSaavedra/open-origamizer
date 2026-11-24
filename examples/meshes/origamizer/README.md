# Origamizer sample models (local only)

This folder holds the sample models distributed with Tomohiro Tachi's
original Origamizer, used here as test inputs for comparison with the paper
(TVCG 2010, Fig. 22 and 23). They were made in Metasequoia, so they are input
meshes, not files produced by Origamizer.

They are **not committed**: they come with Origamizer under its license
(non-commercial use with attribution), and `bunny170` derives from the
Stanford bunny, whose terms forbid commercial use without permission. Both
conflict with redistribution in this GPL repository. Git ignores everything
here except this README.

To use them, copy the files into this folder under these names:

| File | Faces | Shape |
|---|---|---|
| `bell8x8.obj` | 64 quads, 8 x 8 grid | a bell (Gaussian bump); the paper's model (b) |
| `hypar8x8.obj` | 64 quads, 8 x 8 grid | a hyperbolic paraboloid; the paper's model (a) |
| `test_surface_triang.obj` | 512 triangles, 16 x 16 grid | a wavy test surface |
| `bunny170.obj` | 327 triangles | a coarse open bunny (171 vertices, one boundary loop) |

The quads of the two grids are not planar, so load them with
`origamizer.load_mesh(path, triangulate_nonplanar=True)`. When the files are
present, `tests/python/test_origamizer_samples.py` runs on them and
`examples/plot_cut.py` and `examples/plot_layout.py` draw them.
