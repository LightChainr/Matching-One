
from __future__ import annotations
import math
from pathlib import Path
import sys
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p48_sprime_prospective import (  # noqa: E402
    MODEL_ORDER,
    POWER,
    TARGET_SIZES,
    score,
    source_prediction,
    validate_artifact,
)


ARTIFACT_PATH = ROOT / "predictions" / "p48_sprime_correction_20260828.yaml"


def load_artifact() -> dict:
    with ARTIFACT_PATH.open(encoding="utf-8") as handle:
        payload = yaml.safe_load(handle)
    if not isinstance(payload, dict):
        raise ValueError("artifact must be a mapping")
    return payload


class P48ProspectiveSPrimeScorerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.artifact = load_artifact()
        validate_artifact(self.artifact)

    def test_training_excludes_retrospective_holdout(self) -> None:
        training = self.artifact["training_only"]
        self.assertEqual(training["sizes"], [65, 85, 130])
        self.assertEqual(training["excluded_from_fit"], [145, 170])

    def test_model_order_is_frozen(self) -> None:
        self.assertEqual(
            MODEL_ORDER,
            (
                "pure_power_baseline",
                "rank2_log_primary_correction",
                "analytic_inverse_N_competitor",
                "zero_effect",
            ),
        )

    def test_synthetic_log_target_scores_log_at_zero(self) -> None:
        means_y, _ = source_prediction(
            "rank2_log_primary_correction", self.artifact
        )
        scales = [n**POWER for n in TARGET_SIZES]
        raw = [means_y[i] / scales[i] for i in range(2)]
        target = {
            "sizes": list(TARGET_SIZES),
            "independent_of_training_source": True,
            "P4_S_prime": raw,
            "covariance_P4_S_prime": [[1.0e-8, 0.0], [0.0, 1.0e-8]],
        }
        result = score(target, self.artifact)
        rows = {row["name"]: row for row in result["results"]}
        self.assertEqual(result["model_order"], list(MODEL_ORDER))
        self.assertAlmostEqual(
            rows["rank2_log_primary_correction"]["chi_square"],
            0.0,
            places=24,
        )
        self.assertGreater(rows["pure_power_baseline"]["chi_square"], 0.0)
        self.assertGreater(rows["analytic_inverse_N_competitor"]["chi_square"], 0.0)
        self.assertGreater(rows["zero_effect"]["chi_square"], 0.0)

    def test_raw_covariance_is_scaled_to_Y_space(self) -> None:
        target = {
            "sizes": list(TARGET_SIZES),
            "independent_of_training_source": True,
            "P4_S_prime": [0.001, 0.002],
            "covariance_P4_S_prime": [[4.0e-8, 1.0e-8], [1.0e-8, 9.0e-8]],
        }
        result = score(target, self.artifact)
        covariance = result["target_covariance_Y"]
        scales = [n**POWER for n in TARGET_SIZES]
        self.assertAlmostEqual(covariance[0][0], 4.0e-8 * scales[0] ** 2)
        self.assertAlmostEqual(covariance[0][1], 1.0e-8 * scales[0] * scales[1])
        self.assertAlmostEqual(covariance[1][1], 9.0e-8 * scales[1] ** 2)

    def test_source_reuse_and_wrong_sizes_are_rejected(self) -> None:
        base = {
            "sizes": list(TARGET_SIZES),
            "independent_of_training_source": True,
            "P4_S_prime": [0.001, 0.002],
            "covariance_P4_S_prime": [[1.0e-8, 0.0], [0.0, 1.0e-8]],
        }
        reused = dict(base)
        reused["independent_of_training_source"] = False
        with self.assertRaisesRegex(ValueError, "independent"):
            score(reused, self.artifact)

        wrong = dict(base)
        wrong["sizes"] = [170, 265]
        with self.assertRaisesRegex(ValueError, "exactly"):
            score(wrong, self.artifact)

    def test_df2_survival_formula(self) -> None:
        target = {
            "sizes": list(TARGET_SIZES),
            "independent_of_training_source": True,
            "P4_S_prime": [0.001, 0.002],
            "covariance_P4_S_prime": [[1.0e-8, 0.0], [0.0, 1.0e-8]],
        }
        result = score(target, self.artifact)
        for row in result["results"]:
            self.assertEqual(row["df"], 2)
            self.assertAlmostEqual(
                row["chi_square_survival_df2"],
                math.exp(-0.5 * row["chi_square"]),
            )


if __name__ == "__main__":
    unittest.main()
