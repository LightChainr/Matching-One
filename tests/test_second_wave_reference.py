#!/usr/bin/env python3
"""Exact regression contracts for the second-wave matching reference code."""

from __future__ import annotations

from pathlib import Path
import sys
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import exact_matching_polynomial as polynomials  # noqa: E402
import matched_torus_reference as reference  # noqa: E402


class SecondWaveReferenceTests(unittest.TestCase):
    def setUp(self) -> None:
        mp.mp.dps = 90

    def test_finite_matching_identity_is_exact_on_two_topologies(self) -> None:
        p = mp.mpf("0.37")
        geometries = (
            reference.axis_geometry(3),
            reference.diamond_geometry(2),
        )
        for geometry in geometries:
            with self.subTest(geometry=geometry.name):
                result = reference.exact_check(geometry, p)
                self.assertLess(abs(result["difference"]), mp.mpf("1e-70"))

    def test_axis_and_diamond_polynomial_coefficients_are_frozen(self) -> None:
        axis = polynomials.bernstein_to_power(
            polynomials.bernstein_counts(reference.axis_geometry(3))
        )
        diamond = polynomials.bernstein_to_power(
            polynomials.bernstein_counts(reference.diamond_geometry(2))
        )
        self.assertEqual(axis, [-1, 0, 0, 6, 0, 0, 0, -18, 18, -4])
        self.assertEqual(diamond, [-1, 0, 0, 0, 28, -48, 24, 0, -2])

    def test_geometry_contracts_and_exact_size_guard(self) -> None:
        axis = reference.axis_geometry(4)
        diamond = reference.diamond_geometry(3)
        self.assertEqual(axis.n, 16)
        self.assertEqual(len(axis.primal_edges), 2 * axis.n)
        self.assertEqual(len(axis.matching_edges), 4 * axis.n)
        self.assertEqual(diamond.n, 18)
        self.assertEqual(len(diamond.primal_edges), 2 * diamond.n)
        self.assertEqual(len(diamond.matching_edges), 4 * diamond.n)

        with self.assertRaises(ValueError):
            reference.exact_check(reference.axis_geometry(5), mp.mpf("0.37"))


if __name__ == "__main__":
    unittest.main()
