from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_threshold_quantile_certificate import (  # noqa: E402
    build_artifact,
    quantile_bracket,
)


class ExactThresholdQuantileCertificateTests(unittest.TestCase):
    def test_checked_artifact_reproduces(self) -> None:
        artifact = build_artifact()
        checked = json.loads(
            (ROOT / "analysis" / "exact_threshold_quantile_certificate.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(artifact, checked)

    def test_three_frozen_quantiles_are_certified(self) -> None:
        artifact = build_artifact()
        self.assertEqual([row["target"] for row in artifact["quantiles"]], ["1/4", "1/2", "3/4"])
        self.assertTrue(all(row["endpoint_signs_certified"] for row in artifact["quantiles"]))
        self.assertTrue(artifact["quartile_reflection_certified"])
        self.assertEqual(artifact["median_exact"], "1/2")
        self.assertEqual(
            artifact["quantiles"][1],
            {
                "target": "1/2",
                "left": "1/2",
                "right": "1/2",
                "width": "0",
                "left_cdf": "1/2",
                "right_cdf": "1/2",
                "endpoint_signs_certified": True,
            },
        )

    def test_quartile_brackets_reflect_exactly(self) -> None:
        artifact = build_artifact()
        lower, _, upper = artifact["quantiles"]
        lower_left = Fraction(lower["left"])
        lower_right = Fraction(lower["right"])
        upper_left = Fraction(upper["left"])
        upper_right = Fraction(upper["right"])
        self.assertEqual(lower_left + upper_right, 1)
        self.assertEqual(lower_right + upper_left, 1)
        self.assertLessEqual(lower_right - lower_left, Fraction(1, 1 << 24))

    def test_invalid_target_and_nonunique_equations_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "strictly"):
            quantile_bracket([0, 1], Fraction(0))
        with self.assertRaisesRegex(ValueError, "strictly"):
            quantile_bracket([0, 1], Fraction(1))
        # F(x)=4x(1-x) crosses q=1/2 twice on [0,1].
        with self.assertRaisesRegex(ValueError, "exactly one"):
            quantile_bracket([0, 4, -4], Fraction(1, 2), bits=12)
        # Constant CDF candidate has no root.
        with self.assertRaisesRegex(ValueError, "exactly one"):
            quantile_bracket([0], Fraction(1, 2), bits=12)

    def test_linear_control_is_exact(self) -> None:
        self.assertEqual(
            quantile_bracket([0, 1], Fraction(3, 8), bits=12),
            (Fraction(3, 8), Fraction(3, 8)),
        )


if __name__ == "__main__":
    unittest.main()

