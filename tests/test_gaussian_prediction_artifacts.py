from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest

import mpmath as mp
import yaml


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gaussian_semigroup_design import Gaussian, lineage_payload  # noqa: E402


PREDICTIONS = ROOT / "predictions"


def load_yaml(name: str) -> dict:
    with (PREDICTIONS / name).open(encoding="utf-8") as handle:
        payload = yaml.safe_load(handle)
    if not isinstance(payload, dict):
        raise ValueError("prediction artifact must contain a mapping")
    return payload


def fraction_from_mapping(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


class GaussianPredictionArtifactTests(unittest.TestCase):
    def setUp(self) -> None:
        mp.mp.dps = 80

    def test_third_doubling_artifact_matches_exact_generator(self) -> None:
        frozen = load_yaml("gaussian_doubling_third_lineage_20260828.yaml")
        self.assertEqual(frozen["status"], "preregistered_before_target_run")
        self.assertEqual(frozen["claim_level"], "C0")
        generated = lineage_payload(
            Gaussian(*frozen["parent"]["first"]),
            Gaussian(*frozen["parent"]["second"]),
            Gaussian(*frozen["multiplier"]["gaussian_integer"]),
        )
        self.assertEqual(
            generated["parent"]["first"]["translation_group"],
            frozen["parent"]["translation_group"],
        )
        self.assertEqual(
            generated["child"]["first_canonical"]["translation_group"],
            frozen["child"]["translation_group"],
        )
        self.assertEqual(
            generated["child"]["first_canonical"]["pair"],
            frozen["child"]["canonical_D4_representatives"]["first"],
        )
        self.assertEqual(
            generated["child"]["second_canonical"]["pair"],
            frozen["child"]["canonical_D4_representatives"]["second"],
        )
        h4 = generated["harmonic_predictions"]["H4"]
        self.assertEqual(
            Fraction(
                h4["parent_delta"]["numerator"],
                h4["parent_delta"]["denominator"],
            ),
            fraction_from_mapping(frozen["parent"]["exact_delta_cos4"]),
        )
        self.assertEqual(
            Fraction(
                h4["child_delta"]["numerator"],
                h4["child_delta"]["denominator"],
            ),
            fraction_from_mapping(frozen["child"]["exact_delta_cos4_lineage"]),
        )
        target = -mp.power(2, -mp.mpf(13) / 8)
        self.assertLess(
            abs(mp.mpf(str(frozen["primary_prediction"]["target_decimal"])) - target),
            mp.mpf("1e-48"),
        )
        slope_target = mp.power(2, mp.mpf(3) / 8)
        self.assertLess(
            abs(
                mp.mpf(
                    str(
                        frozen["linked_full_curve_predictions"]["mean_Mprime_ratio"][
                            "target_decimal"
                        ]
                    )
                )
                - slope_target
            ),
            mp.mpf("1e-48"),
        )

    def test_norm5_artifact_matches_both_exact_lineages(self) -> None:
        frozen = load_yaml("gaussian_norm5_harmonic_discrimination_20260828.yaml")
        self.assertEqual(frozen["status"], "preregistered_before_child_runs")
        self.assertEqual(frozen["claim_level"], "C0")
        expected_ratios = {
            "H4": Fraction(-14, 25),
            "H8": Fraction(-1054, 625),
            "H12": Fraction(23506, 15625),
        }
        for lineage in frozen["lineages"]:
            generated = lineage_payload(
                Gaussian(*lineage["parent"]["first"]),
                Gaussian(*lineage["parent"]["second"]),
                Gaussian(*lineage["multiplier"]),
            )
            self.assertEqual(
                generated["parent"]["first"]["translation_group"],
                lineage["parent"]["translation_group"],
            )
            self.assertEqual(
                generated["child"]["first_canonical"]["translation_group"],
                lineage["child_canonical_lineage"]["translation_group"],
            )
            self.assertEqual(
                generated["child"]["first_canonical"]["pair"],
                lineage["child_canonical_lineage"]["first"],
            )
            self.assertEqual(
                generated["child"]["second_canonical"]["pair"],
                lineage["child_canonical_lineage"]["second"],
            )
            for harmonic, expected in expected_ratios.items():
                ratio = generated["harmonic_predictions"][harmonic]["angular_ratio"]
                self.assertEqual(
                    Fraction(ratio["numerator"], ratio["denominator"]),
                    expected,
                )

        radial = mp.power(5, -mp.mpf(13) / 8)
        for harmonic, angular in expected_ratios.items():
            target = mp.mpf(angular.numerator) / angular.denominator * radial
            reported = mp.mpf(
                str(frozen["exact_harmonic_predictions"][harmonic]["target_decimal"])
            )
            self.assertLess(abs(reported - target), mp.mpf("1e-29"))


if __name__ == "__main__":
    unittest.main()
