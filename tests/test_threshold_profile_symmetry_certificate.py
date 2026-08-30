from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from threshold_profile_symmetry_certificate import (  # noqa: E402
    build_artifact,
    cdf_has_complement_reflection,
    certify_unique_midpoint_mode,
    density_is_reflection_symmetric,
    reflection_coefficients,
    weights_are_reflection_symmetric,
)


class ThresholdProfileSymmetryCertificateTests(unittest.TestCase):
    def test_checked_artifact_reproduces(self) -> None:
        artifact = build_artifact()
        checked = json.loads(
            (ROOT / "analysis" / "threshold_profile_symmetry_certificate.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(artifact, checked)

    def test_frozen_profile_has_exact_unique_center(self) -> None:
        artifact = build_artifact()
        self.assertTrue(all(artifact["reflection_gates"].values()))
        self.assertTrue(artifact["unique_median_certified"])
        self.assertTrue(artifact["unique_mode_certified"])
        self.assertEqual(artifact["median"], "1/2")
        self.assertEqual(artifact["mode"], "1/2")
        self.assertEqual(artifact["derivative_power_coefficients"], ["3", "-6"])

    def test_reflection_transform_is_exact(self) -> None:
        self.assertEqual(
            reflection_coefficients([Fraction(1), Fraction(2), Fraction(3)]),
            [Fraction(6), Fraction(-8), Fraction(3)],
        )
        self.assertTrue(density_is_reflection_symmetric([Fraction(1, 2), 3, -3, 0]))
        self.assertTrue(cdf_has_complement_reflection([0, Fraction(1, 2), Fraction(3, 2), -1, 0]))

    def test_asymmetric_weights_are_not_certified(self) -> None:
        self.assertTrue(weights_are_reflection_symmetric([1, 2, 2, 1]))
        self.assertFalse(weights_are_reflection_symmetric([1, 2, 3, 1]))

    def test_unique_mode_gate_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "linear derivative"):
            certify_unique_midpoint_mode([Fraction(1)])
        with self.assertRaisesRegex(ValueError, "linear derivative"):
            certify_unique_midpoint_mode([1, 0, -1, 1])
        with self.assertRaisesRegex(ValueError, "positive to negative"):
            certify_unique_midpoint_mode([1, -1, 1])
        with self.assertRaisesRegex(ValueError, "reflection midpoint"):
            certify_unique_midpoint_mode([1, 2, -1])


if __name__ == "__main__":
    unittest.main()

