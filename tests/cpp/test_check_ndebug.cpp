// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

// Release builds define NDEBUG, which removes assert(). This file defines it
// explicitly to pin down that ORIGAMIZER_CHECK still fires without it.
#ifndef NDEBUG
#define NDEBUG 1
#endif

#include <catch2/catch_test_macros.hpp>

#include "origamizer/check.hpp"

TEST_CASE("ORIGAMIZER_CHECK survives NDEBUG") {
  REQUIRE_THROWS_AS(ORIGAMIZER_CHECK(false, "still checked"), origamizer::InvariantError);
}
