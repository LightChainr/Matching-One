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

    def test_norm5_radial_competitors_match_exact_formula(self) -> None:
        frozen = load_yaml("gaussian_norm5_radial_competitors_20260828.yaml")
        self.assertEqual(frozen["claim_level"], "C0")
        angular = -mp.mpf(14) / 25
        expected = {
            "x21_over_4_thermal_level4": Fraction(13, 8),
            "x14_over_3_V13_parity_failure": Fraction(4, 3),
            "x17_over_4_W22_log_leakage": Fraction(9, 8),
        }
        for name, alpha in expected.items():
            model = frozen["models"][name]
            self.assertEqual(Fraction(model["exponent_alpha_in_N"]), alpha)
            target = angular * mp.power(
                5, -mp.mpf(alpha.numerator) / alpha.denominator
            )
            reported = mp.mpf(model["deltaM_child_over_parent"])
            self.assertLess(abs(reported - target), mp.mpf("1e-48"))

        primary = frozen["models"]["x21_over_4_thermal_level4"]
        self.assertEqual(
            primary["full_curve_if_run"]["root_gap_child_over_parent_expression"],
            "-14/625",
        )
        self.assertEqual(
            mp.mpf(str(primary["full_curve_if_run"]["root_gap_child_over_parent"])),
            -mp.mpf(14) / 625,
        )

    def test_cross_norm_residual_transfer_artifact_matches_formula(self) -> None:
        frozen = load_yaml("gaussian_norm2_norm5_residual_transfer_20260828.yaml")
        self.assertEqual(frozen["claim_level"], "C0")
        alpha = mp.mpf(13) / 8
        r2 = -mp.mpf(1)
        r5 = -mp.mpf(14) / 25

        log_expected = (
            r5 * mp.power(5, -alpha) * mp.log(5)
            / (r2 * mp.power(2, -alpha) * mp.log(2))
        )
        log_reported = mp.mpf(
            frozen["transfer_laws"]["logarithmic_Jordan_rank2"]["E5_over_E2"]
        )
        self.assertLess(abs(log_reported - log_expected), mp.mpf("1e-48"))

        for q_length, key in (
            (2, "relative_length_q2"),
            (3, "relative_length_q3"),
            (4, "relative_length_q4"),
            (6, "relative_length_q6"),
        ):
            beta = mp.mpf(q_length) / 2
            expected = (
                r5
                * mp.power(5, -alpha)
                * (mp.power(5, -beta) - 1)
                / (
                    r2
                    * mp.power(2, -alpha)
                    * (mp.power(2, -beta) - 1)
                )
            )
            reported = mp.mpf(frozen["transfer_laws"][key]["E5_over_E2"])
            self.assertLess(abs(reported - expected), mp.mpf("1e-48"))

        for size in ("N65", "N85"):
            source = frozen["source_norm2_residuals"][size]
            e2 = mp.mpf(source["E2"])
            se2 = mp.mpf(source["E2_standard_error"])
            for key, factor in (
                ("logarithmic_Jordan_rank2", log_expected),
                *(
                    (
                        name,
                        mp.mpf(frozen["transfer_laws"][name]["E5_over_E2"]),
                    )
                    for name in (
                        "relative_length_q2",
                        "relative_length_q3",
                        "relative_length_q4",
                        "relative_length_q6",
                    )
                ),
            ):
                target = frozen["transfer_laws"][key]["prospective_E5_from_source"][
                    size
                ]
                self.assertLess(
                    abs(mp.mpf(target["mean"]) - factor * e2),
                    mp.mpf("1e-34"),
                )
                self.assertLess(
                    abs(mp.mpf(target["source_only_standard_error"]) - factor * se2),
                    mp.mpf("1e-34"),
                )


if __name__ == "__main__":
    unittest.main()
