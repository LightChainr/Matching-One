from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from crosswalk_natural_current_geometry import (  # noqa: E402
    crosswalk,
    fit_model,
    gaussian_chi4,
    load_scale,
)


class NaturalCurrentGeometryCrosswalkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = crosswalk(ROOT)

    def test_exact_H4_descriptors(self) -> None:
        self.assertEqual(str(gaussian_chi4(8, 1)[0]), "3713/4225")
        self.assertEqual(str(gaussian_chi4(7, 4)[0]), "-2047/4225")
        self.assertEqual(str(gaussian_chi4(12, 1)[0]), "19873/21025")
        self.assertEqual(str(gaussian_chi4(9, 8)[0]), "-20447/21025")
        for scale in self.payload["data"]["scales"]:
            for descriptor in scale["descriptors"]:
                self.assertTrue(descriptor["z_axis_chi4"]["unit_norm_exact"])
                self.assertEqual(descriptor["tau"], "i")

    def test_full_covariance_and_heldout_fit_contract(self) -> None:
        covariance = self.payload["data"]["covariance"]
        self.assertEqual(len(covariance), 6)
        self.assertTrue(all(len(row) == 6 for row in covariance))
        for model in self.payload["models"]:
            self.assertEqual(model["fit_scales"], [65, 85])
            self.assertEqual(model["heldout_N"], 145)

    def test_H4_geometry_predicts_heldout_rotation(self) -> None:
        models = {row["name"]: row for row in self.payload["models"]}
        pure = models["pure_N_law"]
        h4 = models["one_H4_geometry_covector"]
        extended = models["H4_geometry_plus_A_projective_scalar"]
        self.assertAlmostEqual(h4["training_quadratic"], 2.57523951510026)
        self.assertAlmostEqual(h4["heldout_predictive_quadratic"], 1.5191488422420671)
        self.assertAlmostEqual(h4["heldout_pair_contrast"]["z"], 1.1187917811670134)
        self.assertGreater(pure["training_quadratic"], 23.0)
        self.assertLess(
            h4["training_quadratic"] - extended["training_quadratic"], 0.13
        )
        decomposition = self.payload["diagnosis"]["central_rebound_decomposition"]
        self.assertAlmostEqual(
            decomposition["central_geometry_fraction_of_apparent_rebound"],
            0.22770036221460632,
        )
        self.assertLess(
            decomposition["geometry_aware_N85_anchored"]["quadratic"], 2.12
        )

    def test_N145_values_do_not_select_direction(self) -> None:
        training = [load_scale(ROOT, 65), load_scale(ROOT, 85)]
        heldout = load_scale(ROOT, 145)
        first = fit_model(
            "H4", ["H4_geometry"], training, heldout
        )
        altered = dict(heldout)
        altered["values"] = [100.0, -100.0]
        second = fit_model(
            "H4", ["H4_geometry"], training, altered
        )
        self.assertEqual(first["parameters"], second["parameters"])
        self.assertNotEqual(
            first["heldout_predictive_quadratic"],
            second["heldout_predictive_quadratic"],
        )

    def test_next_geometry_is_exact_angle_flip(self) -> None:
        next_geometry = self.payload["next_same_lineage_geometry"]
        self.assertEqual(next_geometry["child"]["N"], 170)
        parent = [(9, 2), (7, 6)]
        child = [(11, 7), (13, 1)]
        for source, target in zip(parent, child):
            self.assertEqual(gaussian_chi4(*target)[0], -gaussian_chi4(*source)[0])
            self.assertEqual(
                abs(gaussian_chi4(*target)[1]), abs(gaussian_chi4(*source)[1])
            )


if __name__ == "__main__":
    unittest.main()
