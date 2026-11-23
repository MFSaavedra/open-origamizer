# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Matías Saavedra

"""Smoke tests for the compiled core and its exception translation."""

import re

import pytest

import origamizer
from origamizer import _core


def test_cgal_version_is_at_least_5_6():
    match = re.match(r"(\d+)\.(\d+)", origamizer.cgal_version())
    assert match is not None
    assert (int(match.group(1)), int(match.group(2))) >= (5, 6)


def test_failed_check_raises_invariant_error():
    with pytest.raises(origamizer.InvariantError, match="raised from the test") as info:
        _core._fail_check_for_testing("raised from the test")
    assert isinstance(info.value, RuntimeError)
