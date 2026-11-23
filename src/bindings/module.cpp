// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

// Python bindings for the C++ core: the extension module origamizer._core.

#include <CGAL/version.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <string>

#include "origamizer/check.hpp"

namespace nb = nanobind;

NB_MODULE(_core, m) {
  m.doc() = "C++ core of Open Origamizer.";

  nb::exception<origamizer::InvariantError>(m, "InvariantError", PyExc_RuntimeError);

  m.def(
      "cgal_version", [] { return std::string(CGAL_VERSION_STR); },
      "Version of the CGAL headers the core was compiled against.");

  m.def(
      "_fail_check_for_testing",
      [](const std::string& message) { ORIGAMIZER_CHECK(false, message); },
      nb::arg("message"),
      "Fail an invariant check on purpose. Test hook for the InvariantError translation.");
}
