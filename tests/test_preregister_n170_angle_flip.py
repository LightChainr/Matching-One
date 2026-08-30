from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from preregister_n170_angle_flip import preregister  # noqa: E402


class N170AngleFlipPreregistrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = preregister(
            ROOT / "results/p337-natural-current-geometry-crosswalk/latest.json",
            ROOT / "results/server-20260830/P337-natural-current-third-scale-N145/score.json",
        )

    def test_exact_child_geometry(self) -> None:
        geometry = self.payload["geometry"]
        self.assertEqual(geometry["N"], 170)
        self.assertTrue(geometry["exact_angle_flip"])
        self.assertEqual(geometry["H4_covectors_exact"], ["-4633/7225", "6887/7225"])
        self.assertEqual(geometry["charged_projective_scalar"]["value"], 0.5)

    def test_frozen_H4_prediction(self) -> None:
        prediction = self.payload["frozen_H4_only_prediction"]
        self.assertAlmostEqual(prediction["absolute_K_A"][0], 0.0039355336023874695)
        self.assertAlmostEqual(prediction["absolute_K_A"][1], -0.0058502093502358085)
        self.assertAlmostEqual(prediction["pair_second_minus_first"], -0.009785742952623279)
        self.assertLess(prediction["H4_amplitude"], 0.0)
        self.assertEqual(prediction["A_projective_scalar"], 0.0)

    def test_power_and_production_freeze(self) -> None:
        production = self.payload["production"]
        self.assertEqual(production["samples_per_shape"], 8_000_000)
        self.assertEqual(production["batches"], 80)
        self.assertEqual(production["seed"], 202608337170)
        self.assertEqual(production["machine"]["id"], "033945d8bf8b47a7acf475c595169e07")
        self.assertGreater(self.payload["power"]["projected_pair_z_at_H4_target"], 5.6)
        self.assertLess(
            self.payload["power"]["required_samples_for_5sigma_pair_vs_zero"],
            production["samples_per_shape"],
        )

    def test_no_harmonic_revoting(self) -> None:
        contract = self.payload["scoring_contract"]
        self.assertTrue(contract["no_H4_H8_vote"])
        self.assertTrue(contract["no_exponent_fit"])
        self.assertTrue(contract["no_post_reveal_basis_change"])


if __name__ == "__main__":
    unittest.main()
