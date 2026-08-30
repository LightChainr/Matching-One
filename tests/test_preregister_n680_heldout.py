from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from preregister_n680_heldout import preregister  # noqa: E402


class N680HeldoutPreregistrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = preregister(
            ROOT / "results/p337-three-generation-recurrence/latest.json",
            ROOT / "results/server-20260830/P337-N340-second-child/score.json",
        )

    def test_exact_child(self) -> None:
        geometry = self.payload["geometry"]
        self.assertEqual(geometry["N"], 680)
        self.assertEqual(geometry["H4_covectors_exact"], ["-4633/7225", "6887/7225"])
        self.assertEqual(geometry["Smith_classes"], [[2, 340], [2, 340]])
        self.assertTrue(geometry["exact_angle_flip"])

    def test_four_predictions_are_frozen(self) -> None:
        models = self.payload["frozen_models_in_scoring_order"]
        self.assertEqual([model["name"] for model in models], [
            "two_mode_recurrence", "single_frozen_lambda0", "single_free_lambda", "scale_neutral",
        ])
        self.assertAlmostEqual(models[0]["H4_amplitude"], -0.001841297042231899)
        self.assertAlmostEqual(models[2]["H4_amplitude"], -0.003037819395650872)
        self.assertTrue(all(model["pair_second_minus_first"] < 0.0 for model in models))

    def test_power_contract(self) -> None:
        production = self.payload["production"]
        self.assertEqual(production["samples_per_shape"], 120_000_000)
        self.assertEqual(production["batches"], 80)
        self.assertEqual(production["seed"], 202608337680)
        power = self.payload["power"]
        self.assertGreater(power["projected_measurement_z_for_gaps"]["two_mode_vs_free_single"], 3.0)
        self.assertGreater(power["projected_measurement_z_for_gaps"]["two_mode_vs_scale_neutral"], 17.0)
        self.assertLess(power["required_samples_for_3sigma_two_mode_vs_free_single"], 120_000_000)

    def test_no_refit(self) -> None:
        contract = self.payload["scoring_contract"]
        self.assertTrue(contract["projective_scalar_zero_control"])
        self.assertTrue(contract["no_model_refit"])
        self.assertTrue(contract["no_exponent_fit"])
        self.assertTrue(contract["no_post_reveal_basis_change"])


if __name__ == "__main__":
    unittest.main()
