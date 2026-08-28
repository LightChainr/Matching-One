from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p48_sprime_frozen import (  # noqa: E402
    basis,
    prediction_covariance,
    score,
    validate_manifest,
)


PREDICTIONS = ROOT / "predictions"


def load_yaml(name: str) -> dict:
    with (PREDICTIONS / name).open(encoding="utf-8") as handle:
        payload = yaml.safe_load(handle)
    if not isinstance(payload, dict):
        raise ValueError("prediction artifact must contain a mapping")
    return payload


class P48FrozenSPrimeScorerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.manifest = load_yaml("p48_sprime_scoring_manifest_20260828.yaml")
        self.sizes, self.power, self.models = validate_manifest(self.manifest)
        self.by_name = {model["name"]: model for model in self.models}

    def test_manifest_matches_separate_correction_artifacts(self) -> None:
        q2 = load_yaml("p48_sprime_q2_correction_20260828.yaml")
        log = load_yaml("p48_sprime_jordan_log_20260828.yaml")
        q2_model = self.by_name["q2_even_scalar"]
        log_model = self.by_name["rank2_jordan_log"]

        self.assertEqual(
            q2_model["parameters"],
            [q2["fit"]["A"], q2["fit"]["B"]],
        )
        self.assertEqual(
            q2_model["parameter_covariance"],
            q2["fit"]["covariance_AB"],
        )
        self.assertEqual(
            log_model["parameters"],
            [log["fit"]["A"], log["fit"]["B_logN"]],
        )
        self.assertEqual(
            log_model["parameter_covariance"],
            log["fit"]["covariance_AB"],
        )

    def test_frozen_predictions_are_recomputed_from_parameters(self) -> None:
        for model in self.models:
            frozen = {
                int(row["N"]): row for row in model["frozen_predictions"]
            }
            if model["basis"] == "zero":
                for n in self.sizes:
                    self.assertEqual(frozen[n]["P4_S_prime"], 0.0)
                    self.assertEqual(frozen[n]["source_fit_se"], 0.0)
                continue

            design = [basis(n, model, self.power) for n in self.sizes]
            parameters = [float(value) for value in model["parameters"]]
            covariance = prediction_covariance(
                design, model["parameter_covariance"]
            )
            for index, n in enumerate(self.sizes):
                predicted = math.fsum(
                    coefficient * parameter
                    for coefficient, parameter in zip(
                        design[index], parameters
                    )
                )
                self.assertAlmostEqual(
                    predicted,
                    float(frozen[n]["P4_S_prime"]),
                    places=15,
                )
                self.assertAlmostEqual(
                    math.sqrt(covariance[index][index]),
                    float(frozen[n]["source_fit_se"]),
                    places=15,
                )

    def test_synthetic_q2_target_scores_q2_at_zero_without_refit(self) -> None:
        frozen = self.by_name["q2_even_scalar"]["frozen_predictions"]
        target = {
            "sizes": list(self.sizes),
            "independent_of_retrospective_source": True,
            "P4_S_prime": [row["P4_S_prime"] for row in frozen],
            "covariance_P4_S_prime": [
                [1.0e-8, 0.0],
                [0.0, 1.0e-8],
            ],
        }
        result = score(target, self.manifest)
        self.assertEqual(
            result["scoring_order"],
            [
                "pure_N^-5/4",
                "zero_effect",
                "q2_even_scalar",
                "rank2_jordan_log",
            ],
        )
        scores = {row["name"]: row["chi_square"] for row in result["results"]}
        self.assertAlmostEqual(scores["q2_even_scalar"], 0.0, places=24)
        self.assertGreater(scores["pure_N^-5/4"], 0.0)
        self.assertGreater(scores["zero_effect"], 0.0)
        self.assertGreater(scores["rank2_jordan_log"], 0.0)

    def test_source_reuse_is_rejected(self) -> None:
        target = {
            "sizes": list(self.sizes),
            "independent_of_retrospective_source": False,
            "P4_S_prime": [0.0, 0.0],
            "covariance_P4_S_prime": [
                [1.0e-8, 0.0],
                [0.0, 1.0e-8],
            ],
        }
        with self.assertRaisesRegex(ValueError, "independent"):
            score(target, self.manifest)


if __name__ == "__main__":
    unittest.main()
