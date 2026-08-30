from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from preregister_n340_second_child import preregister  # noqa: E402


class N340SecondChildPreregistrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = preregister(
            ROOT / "results/server-20260830/P337-natural-current-scale-N85/score.json",
            ROOT / "results/server-20260830/P337-N170-angle-flip/score.json",
        )

    def test_exact_second_child_geometry(self) -> None:
        geometry = self.payload["geometry"]
        self.assertEqual(geometry["N"], 340)
        self.assertTrue(geometry["exact_angle_flip"])
        self.assertEqual(geometry["H4_covectors_exact"], ["4633/7225", "-6887/7225"])
        self.assertEqual(geometry["Smith_classes"], [[2, 170], [2, 170]])

    def test_frozen_model_order_and_sign(self) -> None:
        models = self.payload["frozen_models_in_scoring_order"]
        self.assertEqual([model["name"] for model in models], [
            "nominal_area_H4", "observed_N85_to_N170_effective", "scale_neutral",
        ])
        self.assertTrue(all(model["pair_second_minus_first"] > 0 for model in models))
        self.assertAlmostEqual(models[0]["H4_amplitude"], -0.0036024562797008176)
        self.assertAlmostEqual(models[1]["H4_amplitude"], -0.007692100036533405)
        self.assertAlmostEqual(models[2]["H4_amplitude"], -0.011111494145226835)

    def test_power_and_zero_control(self) -> None:
        production = self.payload["production"]
        self.assertEqual(production["samples_per_shape"], 12_000_000)
        self.assertEqual(production["batches"], 80)
        self.assertEqual(production["seed"], 202608337340)
        power = self.payload["power"]["projected_measurement_z_for_gaps"]
        self.assertGreater(power["nominal_vs_effective"], 5.0)
        self.assertGreater(power["effective_vs_neutral"], 4.0)
        self.assertTrue(self.payload["scoring_contract"]["projective_scalar_zero_control"])

    def test_no_post_reveal_fit(self) -> None:
        contract = self.payload["scoring_contract"]
        self.assertTrue(contract["no_exponent_fit"])
        self.assertTrue(contract["no_post_reveal_model_or_basis_change"])
        self.assertTrue(contract["no_H4_H8_vote"])


if __name__ == "__main__":
    unittest.main()
