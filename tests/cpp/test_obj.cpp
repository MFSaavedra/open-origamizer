// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "origamizer/check.hpp"
#include "origamizer/mesh.hpp"
#include "support.hpp"

using Catch::Matchers::ContainsSubstring;
using origamizer::InvariantError;
using origamizer::parse_obj;

namespace {
origamizer::PolygonSoup parse(const std::string& text) {
  std::istringstream in(text);
  return parse_obj(in, "test.obj");
}
}  // namespace

TEST_CASE("parse_obj keeps polygon faces") {
  const auto soup = parse(
      "# a pentagon\nv 0 0 0\nv 1 0 0\nv 1 1 0\nv 0.5 1.5 0\nv 0 1 0\n"
      "f 1 2 3 4 5\n");
  REQUIRE(soup.points.size() == 5);
  REQUIRE(soup.faces == std::vector<std::vector<std::size_t>>{{0, 1, 2, 3, 4}});
  REQUIRE(soup.points[3] == origamizer::Point3(0.5, 1.5, 0));
}

TEST_CASE("parse_obj reads texture and normal references and negative indices") {
  const auto soup = parse("v 0 0 0\nv 1 0 0\nv 0 1 0\nvt 0 0\nvn 0 0 1\nf 1/1/1 2//1 -1/1\n");
  REQUIRE(soup.faces == std::vector<std::vector<std::size_t>>{{0, 1, 2}});
}

TEST_CASE("parse_obj ignores names, groups, materials, smoothing, and CRLF line ends") {
  const auto soup = parse(
      "mtllib m.mtl\r\no thing\r\ng part\r\nusemtl red\r\ns off\r\nv 0 0 0\r\nv 1 0 0 1\r\n"
      "v 0 1 0\r\nf 1 2 3   # trailing comment\r\n");
  REQUIRE(soup.points.size() == 3);
  REQUIRE(soup.faces.size() == 1);
}

TEST_CASE("parse_obj rejects unsupported elements and names the line") {
  REQUIRE_THROWS_WITH(parse("v 0 0 0\nl 1 2\n"),
                      ContainsSubstring("test.obj:2") && ContainsSubstring("'l'"));
}

TEST_CASE("parse_obj rejects bad vertex references") {
  REQUIRE_THROWS_WITH(parse("v 0 0 0\nv 1 0 0\nv 0 1 0\nf 0 1 2\n"),
                      ContainsSubstring("bad vertex reference '0'"));
  REQUIRE_THROWS_WITH(parse("v 0 0 0\nf a b c\n"), ContainsSubstring("bad vertex reference"));
  REQUIRE_THROWS_WITH(parse("v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 9\n"),
                      ContainsSubstring("test.obj:4") && ContainsSubstring("only 3"));
  REQUIRE_THROWS_WITH(parse("v 0 0 0\nf -5 1 1\n"), ContainsSubstring("before the first"));
}

TEST_CASE("parse_obj rejects malformed coordinates") {
  REQUIRE_THROWS_WITH(parse("v 1 2\n"), ContainsSubstring("three coordinates"));
  REQUIRE_THROWS_WITH(parse("v 1 2 x\n"), ContainsSubstring("three coordinates"));
  REQUIRE_THROWS_WITH(parse("v 1 2 inf\n"), ContainsSubstring("not finite"));
  REQUIRE_THROWS_WITH(parse("v 1 2 3 junk\n"), ContainsSubstring("unexpected 'junk'"));
  REQUIRE_THROWS_AS(parse("v 1 2\n"), InvariantError);
}

TEST_CASE("read_obj reads an example mesh and reports a missing file") {
  const auto cube = origamizer::read_obj(test::example_path("cube"));
  REQUIRE(cube.points.size() == 8);
  REQUIRE(cube.faces.size() == 6);
  for (const auto& face : cube.faces) REQUIRE(face.size() == 4);
  REQUIRE_THROWS_AS(origamizer::read_obj("/nonexistent/mesh.obj"), std::runtime_error);
}
