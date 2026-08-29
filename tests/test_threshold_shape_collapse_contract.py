from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import threshold_shape_collapse_contract as contract  # noqa: E402


class ThresholdShapeCollapseContractTests(unittest.TestCase):
    def setUp(self) -> None:
        self.reference = contract.fixture(range(7), (1, 2, 4, 6, 4, 2, 1))

    def test_frozen_quantile_grid(self) -> None:
        self.assertEqual(contract.QUANTILE_GRID, tuple(Fraction(v, 100) for v in (5, 10, 25, 50, 75, 90, 95)))

    def test_positive_affine_transform_has_zero_shape_residual(self) -> None:
        transformed = self.reference.affine(Fraction(5), Fraction(3))
        result = contract.compare_shapes(self.reference, transformed)
        self.assertEqual(result["location_shift"], 11)
        self.assertEqual(result["scale_ratio"], 3)
        self.assertEqual(result["shape_sse"], 0)
        self.assertEqual(result["shape_max_abs"], 0)

    def test_tail_deformation_is_not_absorbed_by_location_and_scale(self) -> None:
        deformed = contract.fixture(range(7), (3, 1, 2, 6, 4, 2, 2))
        result = contract.compare_shapes(self.reference, deformed)
        self.assertEqual(result["location_shift"], 0)
        self.assertEqual(result["scale_ratio"], 1)
        self.assertGreater(result["shape_sse"], 0)
        self.assertGreater(result["shape_max_abs"], 0)

    def test_zero_iqr_fails_closed(self) -> None:
        concentrated = contract.fixture((0, 1, 2), (1, 10, 1))
        with self.assertRaisesRegex(ValueError, "interquartile scale"):
            contract.standardized_profile(concentrated)

    def test_invalid_distribution_and_affine_scale_fail(self) -> None:
        with self.assertRaises(ValueError):
            contract.WeightedDistribution((Fraction(0), Fraction(0)), (1, 1))
        with self.assertRaises(ValueError):
            self.reference.affine(Fraction(0), Fraction(0))

    def test_checked_in_artifacts_reproduce(self) -> None:
        artifact = contract.build_artifact()
        checked_json = json.loads((ROOT / "results/threshold-shape-collapse/latest.json").read_text())
        checked_md = (ROOT / "results/threshold-shape-collapse/latest.md").read_text()
        self.assertEqual(checked_json, artifact)
        self.assertEqual(checked_md, contract.render_markdown(artifact))


if __name__ == "__main__":
    unittest.main()
