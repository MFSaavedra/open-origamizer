# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Matías Saavedra

"""Open Origamizer: fold one uncut sheet into a polyhedral surface by tucking.

The algorithms live in the C++ extension module ``origamizer._core``. This
package adds the file writers, the fabrication helpers, and the command line.
"""

from ._core import InvariantError, cgal_version

__all__ = ["InvariantError", "cgal_version"]
