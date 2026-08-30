#!/usr/bin/env python3
"""Tests for the exact Issue 16 odd-jet invariant primitive."""

from __future__ import annotations

import copy
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import odd_jet_invariants as odd_jet  # noqa: E402


CONTRACT_PATH = ROOT / "analysis" / "odd_jet_invariant_contract.json"


class OddJetInvariantTests(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
        self.orders = (1, 3, 5)
        self.derivatives = (Fraction(2), Fraction(-16), Fraction(384))
        self.covariance = (
            (Fraction(1, 100), Fraction(0), Fraction(0)),
            (Fraction(0), Fraction(4, 25), Fraction(0)),
            (Fraction(0), Fraction(0), Fraction(64)),
        )

    def test_checked_in_contract_is_exact_and_data_free(self) -> None:
        result = odd_jet.validate_contract(self.contract)
        self.assertEqual(result["invariants"], {"3": "-2", "5": "12"})
        self.assertEqual(result["invariant_covariance"], [["37/400", "-9/10"], ["-9/10", "145/16"]])
        self.assertEqual(result["coordinate_rescales_verified"], ["2", "-3/2"])
        self.assertFalse(result["contains_estimate"])
        self.assertEqual(result["parent_issue"], "remain open")

    def test_coordinate_rescaling_preserves_invariants_and_covariance(self) -> None:
        expected_invariants = odd_jet.normalized_invariants(self.orders, self.derivatives)
        expected_covariance = odd_jet.propagate_covariance(self.orders, self.derivatives, self.covariance)
        for scale in (Fraction(7, 3), Fraction(-5, 2)):
            with self.subTest(scale=scale):
                derivatives = odd_jet.rescale_coordinate(self.orders, self.derivatives, scale)
                covariance = odd_jet.rescale_covariance(self.orders, self.covariance, scale)
                self.assertEqual(odd_jet.normalized_invariants(self.orders, derivatives), expected_invariants)
                self.assertEqual(odd_jet.propagate_covariance(self.orders, derivatives, covariance), expected_covariance)

    def test_observable_rescaling_is_not_mislabeled_coordinate_invariance(self) -> None:
        original = odd_jet.normalized_invariants(self.orders, self.derivatives)
        doubled = odd_jet.normalized_invariants(
            self.orders,
            odd_jet.rescale_observable(self.derivatives, Fraction(2)),
        )
        self.assertEqual(doubled[3], original[3] / 4)
        self.assertEqual(doubled[5], original[5] / 16)

    def test_invalid_jet_and_covariance_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "positive odd"):
            odd_jet.normalized_invariants((1, 2, 3), (Fraction(1), Fraction(2), Fraction(3)))
        with self.assertRaisesRegex(ValueError, "first derivative"):
            odd_jet.normalized_invariants((1, 3), (Fraction(0), Fraction(1)))
        changed = [list(row) for row in self.covariance]
        changed[0][1] = Fraction(1)
        with self.assertRaisesRegex(ValueError, "symmetric"):
            odd_jet.propagate_covariance(self.orders, self.derivatives, changed)

    def test_contract_drift_and_empirical_fields_fail(self) -> None:
        changed = copy.deepcopy(self.contract)
        changed["expected"]["invariants"]["3"] = "-5/3"
        with self.assertRaisesRegex(ValueError, "stored normalized invariants drift"):
            odd_jet.validate_contract(changed)

        changed = copy.deepcopy(self.contract)
        changed["synthetic_fixture"]["samples"] = 1000
        with self.assertRaisesRegex(ValueError, "empirical fields"):
            odd_jet.validate_contract(changed)


if __name__ == "__main__":
    unittest.main()
