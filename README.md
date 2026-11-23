# Open Origamizer

An independent, open-source implementation of the Origamizer method: given a polyhedral surface, compute a crease pattern that folds one uncut sheet into it, with the excess material tucked away on one side. The first fabrication target is rigid face tiles 3D printed onto fabric.

**Status:** planning. No code yet.

## Background and clean-room policy

Origamizer is Tomohiro Tachi's method and software. The software is closed source (a Win32 binary, v0.47, 2008), so this project implements the method from the published papers only. Nothing here comes from decompiling or otherwise reverse-engineering the original program. This project is not affiliated with or endorsed by Tomohiro Tachi or the University of Tokyo.

- Tachi. Origamizing Polyhedral Surfaces. IEEE TVCG 16(2):298-311, 2010. [Author PDF](https://origami.c.u-tokyo.ac.jp/~tachi/cg/OrigamizingPolyhedralSurfaces_tachi_IEEETVCG.pdf)
- Demaine, Tachi. Origamizer: A Practical Algorithm for Folding Any Polyhedron. SoCG 2017. [doi:10.4230/LIPIcs.SoCG.2017.34](https://doi.org/10.4230/LIPIcs.SoCG.2017.34) (CC-BY 3.0; a copy is in [docs/papers](docs/papers/))
- Tachi's software: https://tsg.ne.jp/TT/software/

## Plan

1. **Core:** the practical heuristic of TVCG 2010 (tucking molecules, two-step mapping, edge splitting), with [FOLD](https://github.com/edemaine/fold) and SVG export for [Origami Simulator](https://github.com/amandaghassaei/OrigamiSimulator).
2. **Fabrication:** printable face tiles on fabric, beveled so the fabric tucks fit between neighbouring tiles.
3. **Theory (stretch):** the provably correct SoCG 2017 construction.

The implementation will be a C++ core (using CGAL) with Python bindings. Tooling is not set up yet; [CLAUDE.md](CLAUDE.md) holds the working plan.

## License

Copyright (C) 2026 Matías Saavedra

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License in [LICENSE](LICENSE) for more details.

Papers in `docs/papers/` keep their own licenses; see [docs/papers/README.md](docs/papers/README.md).
