#!/usr/bin/env python3
"""Tests for the independent Issue 14 exact interval certificate verifier."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import verify_gadget_comparison_certificate as verifier  # noqa: E402


CERTIFICATE_PATH = ROOT / "analysis" / "synthetic_gadget_comparison_certificate.json"


class GadgetComparisonCertificateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.certificate = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))

    def verify(self, certificate=None):
        return verifier.verify_certificate(self.certificate if certificate is None else certificate)

    def test_checked_in_fixture_verifies_without_claiming_a_bound(self) -> None:
        result = self.verify()
        self.assertEqual(result["status"], "valid_synthetic_fixture_no_bound")
        self.assertEqual([claim["separation"] for claim in result["claims"]], ["1/60", "1/15"])
        self.assertFalse(result["floating_point_used"])
        self.assertFalse(result["proves_new_bound"])
        self.assertEqual(result["parent_issue"], "remain open")

    def test_perturbed_interval_and_strict_boundary_fail(self) -> None:
        changed = copy.deepcopy(self.certificate)
        changed["claims"][0]["lhs"]["lower"] = "7/12"
        with self.assertRaisesRegex(ValueError, "not certified"):
            self.verify(changed)

        touching = copy.deepcopy(self.certificate)
        touching["claims"][0]["lhs"]["lower"] = "7/12"
        touching["claims"][0]["strict"] = False
        result = self.verify(touching)
        self.assertEqual(result["claims"][0]["separation"], "0")

    def test_noncanonical_rationals_and_floats_fail(self) -> None:
        changed = copy.deepcopy(self.certificate)
        changed["claims"][0]["lhs"]["lower"] = "6/10"
        with self.assertRaisesRegex(ValueError, "canonically encoded"):
            self.verify(changed)

        changed = copy.deepcopy(self.certificate)
        changed["claims"][0]["lhs"]["lower"] = 0.6
        with self.assertRaisesRegex(ValueError, "floating point"):
            self.verify(changed)

    def test_theorem_or_new_bound_flags_fail(self) -> None:
        for key in ("theorem_claim", "new_bound_claim"):
            with self.subTest(key=key):
                changed = copy.deepcopy(self.certificate)
                changed["claim_boundary"][key] = True
                with self.assertRaises(ValueError):
                    self.verify(changed)


if __name__ == "__main__":
    unittest.main()
