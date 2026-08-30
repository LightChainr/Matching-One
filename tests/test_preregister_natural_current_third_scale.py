from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from preregister_natural_current_third_scale import preregister  # noqa: E402


class NaturalCurrentThirdScalePreregistrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = preregister(
            ROOT / "analysis/p337_natural_current_scale_preregistration.json",
            ROOT / "results/server-20260830/P337-natural-current-scale-N85/score.json",
        )

    def test_N145_geometry_is_primitive_equal_area(self) -> None:
        geometry = self.payload["geometry"]
        self.assertEqual(geometry["N"], 145)
        self.assertTrue(all(geometry["audit"].values()))
        self.assertEqual(geometry["first"], [12, 1])
        self.assertEqual(geometry["second"], [9, 8])

    def test_four_targets_are_frozen_without_N145(self) -> None:
        targets = self.payload["frozen_targets_at_N145"]
        self.assertAlmostEqual(targets["source_fitted_scale_neutral"]["value"], 0.025592612883702993)
        self.assertAlmostEqual(targets["source_fitted_project_H4"]["value"], 0.010744776300914267)
        self.assertAlmostEqual(targets["secondary_post_reveal_effective_transfer"]["value"], 0.003504048135470623)
        self.assertEqual(targets["secondary_post_reveal_effective_transfer"]["tier"], "secondary")
        self.assertIn("no N145 values", self.payload["claim_boundary"])

    def test_sample_size_covers_primary_curvature_discriminator(self) -> None:
        design = self.payload["design"]
        self.assertEqual(design["samples_per_shape"], 2400000)
        self.assertGreater(
            design["samples_per_shape"],
            design["required_samples_3sigma"]["H4_vs_effective"],
        )


if __name__ == "__main__":
    unittest.main()
