// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Matías Saavedra

// Private helper: disjoint sets (union-find) over 0..n-1.

#pragma once

#include <cstddef>
#include <numeric>
#include <utility>
#include <vector>

namespace origamizer::detail {

/// Union-find with path halving. The smaller index becomes the root of a
/// union, so representatives do not depend on the order of the unions.
class DisjointSets {
 public:
  explicit DisjointSets(std::size_t n) : parent_(n) {
    std::iota(parent_.begin(), parent_.end(), std::size_t{0});
  }

  std::size_t find(std::size_t x) {
    while (parent_[x] != x) {
      parent_[x] = parent_[parent_[x]];
      x = parent_[x];
    }
    return x;
  }

  /// Joins the sets of a and b; returns false if they were already joined.
  bool unite(std::size_t a, std::size_t b) {
    a = find(a);
    b = find(b);
    if (a == b) return false;
    if (b < a) std::swap(a, b);
    parent_[b] = a;
    return true;
  }

 private:
  std::vector<std::size_t> parent_;
};

}  // namespace origamizer::detail
