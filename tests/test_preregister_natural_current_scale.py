from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from preregister_natural_current_scale import preregister  # noqa: E402


class NaturalCurrentScalePreregistrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = preregister(
            ROOT / "results/p337-n65-charged-activity-net/latest.json"
        )

    def test_n65_natural_coordinate(self) -> None:
        natural = self.payload["natural_coordinate"]
        self.assertAlmostEqual(natural["value"][0], -0.04156368093109216)
        self.assertAlmostEqual(natural["value"][1], 0.027916978729333172)
        self.assertAlmostEqual(natural["value"][2], 0.06948065966042533)
        self.assertAlmostEqual(natural["standard_error"][2], 0.020366733051848874)

    def test_three_targets_are_frozen(self) -> None:
        targets = self.payload["frozen_targets_at_N85"]
        self.assertEqual(list(self.payload["scoring_contract"]["model_order"]),
                         ["zero", "source_fitted_scale_neutral", "source_fitted_project_H4"])
        self.assertEqual(targets["zero"]["value"], 0.0)
        self.assertAlmostEqual(targets["source_fitted_scale_neutral"]["value"], 0.06948065966042533)
        self.assertAlmostEqual(targets["source_fitted_project_H4"]["scale_ratio"], 0.6466636515105115)

    def test_N85_design_does_not_use_N85_values(self) -> None:
        design = self.payload["design"]
        self.assertEqual(design["geometry"]["N"], 85)
        self.assertEqual(design["samples_per_shape"], 200000)
        self.assertGreater(
            design["samples_per_shape"],
            design["power_extrapolation"]["required_samples_3sigma"]["scale_neutral_vs_H4"],
        )
        self.assertIn("N85 values absent", self.payload["claim_boundary"])


if __name__ == "__main__":
    unittest.main()
