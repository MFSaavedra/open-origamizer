// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

#include "origamizer/check.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <string>

using Catch::Matchers::ContainsSubstring;

TEST_CASE("ORIGAMIZER_CHECK passes when the condition holds") {
  REQUIRE_NOTHROW(ORIGAMIZER_CHECK(1 + 1 == 2, "arithmetic"));
}

TEST_CASE("ORIGAMIZER_CHECK throws InvariantError with message, condition, and location") {
  REQUIRE_THROWS_AS(ORIGAMIZER_CHECK(1 + 1 == 3, "arithmetic is broken"),
                    origamizer::InvariantError);
  REQUIRE_THROWS_WITH(ORIGAMIZER_CHECK(1 + 1 == 3, "arithmetic is broken"),
                      ContainsSubstring("arithmetic is broken") &&
                          ContainsSubstring("1 + 1 == 3") &&
                          ContainsSubstring("test_check.cpp"));
}

TEST_CASE("ORIGAMIZER_CHECK evaluates its condition exactly once") {
  int calls = 0;
  ORIGAMIZER_CHECK(++calls == 1, "side effect");
  REQUIRE(calls == 1);
}

TEST_CASE("ORIGAMIZER_CHECK accepts a runtime message") {
  const std::string message = "face " + std::to_string(7) + " is not convex";
  REQUIRE_THROWS_WITH(ORIGAMIZER_CHECK(false, message), ContainsSubstring("face 7 is not convex"));
}
