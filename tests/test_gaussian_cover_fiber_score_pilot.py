#!/usr/bin/env python3

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gaussian_cover_fiber_score_pilot import (  # noqa: E402
    LINEAGES,
    Lineage,
    cyclic_label,
    exact_tiny_oracle,
    fiber_parent_label,
    gaussian_integer_torus,
    run_pilot,
    score_values,
    validate_lineage,
)


def test_frozen_norm_two_lineages_have_exact_cyclic_fibers() -> None:
    expected = [29, 24]
    for lineage, multiplier in zip(LINEAGES, expected):
        validate_lineage(lineage)
        assert lineage.parent_multiplier == multiplier
        parent = gaussian_integer_torus(*lineage.parent)
        child = gaussian_integer_torus(*lineage.child)
        assert child.n == 2 * parent.n
        for point in child.coordinates:
            child_label = cyclic_label(*lineage.child, point)
            parent_label = cyclic_label(*lineage.parent, point)
            assert fiber_parent_label(lineage, child_label) == parent_label
            assert fiber_parent_label(lineage, child_label + parent.n) == parent_label


def test_tiny_lineage_and_score_modes() -> None:
    lineage = Lineage((2, 1), (1, 3), 5, 2)
    validate_lineage(lineage)
    labels = [True] * 5 + [False] * 5
    trivial, detail = score_values(labels, 0.5)
    assert trivial == 0.0
    assert detail == 5


def test_exact_score_identity_matches_polynomial_finite_difference() -> None:
    oracle = exact_tiny_oracle()
    assert oracle["p"] == str(Fraction(2, 5))
    assert len(oracle["checks"]) == 4
    assert all(check["equal"] for check in oracle["checks"].values())


def test_pilot_is_worker_invariant_and_has_full_covariance() -> None:
    rows_one, summary_one = run_pilot(80, 4, 0.5, 226, 1)
    rows_two, summary_two = run_pilot(80, 4, 0.5, 226, 2)
    assert rows_one == rows_two
    assert summary_one == summary_two
    assert len(summary_one["primitive_covariance_of_mean"]) == 4
    assert all(len(row) == 4 for row in summary_one["primitive_covariance_of_mean"])
    assert len(summary_one["derived_covariance_of_mean"]) == 4
    assert all(len(row) == 4 for row in summary_one["derived_covariance_of_mean"])
