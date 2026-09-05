
from __future__ import annotations
import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_issue43_secondary import (  # noqa: E402
    P48_ARTIFACT_SHA256,
    X17_ARTIFACT_SHA256,
    score,
    sha256,
)


def primary_score() -> dict:
    return {
        "protocol": "Issue #43 prospective N=185/265 two-spin4 full-curve score",
        "status": "frozen primary score; no refit",
        "prediction_artifact_sha256": (
            "a370e79a10854341fac3ee75e8c518dbf3533e8c077cba2c2ec1018178144f44"
        ),
        "scores": {
            "DeltaM": {
                "observed": [0.00020, 0.00016],
                "sampling_se": [0.00001, 0.000012],
                "target_chi_square": 0.21,
                "target_df": 2,
                "zero_chi_square": 522.0,
                "zero_df": 2,
            },
            "DeltaS": {
                "observed": [0.00007, 0.000069],
                "sampling_se": [0.000008, 0.000009],
                "target_chi_square": 0.08,
                "target_df": 2,
                "zero_chi_square": 106.0,
                "zero_df": 2,
            },
        },
    }


class Issue43SecondaryScoreTests(unittest.TestCase):
    @property
    def x17_artifact(self) -> Path:
        return ROOT / "predictions/x17_spin4_competitor_20260828.yaml"

    @property
    def p48_artifact(self) -> Path:
        return ROOT / "predictions/p48_sprime_correction_20260828.yaml"

    def test_artifacts_are_immutable_inputs(self) -> None:
        self.assertEqual(sha256(self.x17_artifact), X17_ARTIFACT_SHA256)
        self.assertEqual(sha256(self.p48_artifact), P48_ARTIFACT_SHA256)

    def test_fixed_order_exclusion_and_not_scorable_boundary(self) -> None:
        result = score(primary_score(), self.x17_artifact, self.p48_artifact)
        self.assertEqual(
            result["stage_order"],
            [
                "original_x21_H4_two_sector",
                "x17_over_4_H4_adversarial_radial",
                "zero_effect",
                "predeclared_shared_H4_plus_H12",
                "issue72_P48_S_prime",
            ],
        )
        self.assertEqual(result["stages"][1]["status"], "SCORED_FROZEN_NO_REFIT")
        self.assertEqual(result["stages"][1]["target_refit_parameters"], 0)
        self.assertEqual(result["stages"][3]["status"], "NOT_SCORABLE")
        self.assertIn("forbidden target refit", result["stages"][3]["reason"])
        self.assertEqual(result["stages"][4]["status"], "READY_AWAITING_DERIVATIVE_TARGET")
        self.assertEqual(result["excluded_models"][0]["scored"], False)

    def test_p48_is_scored_only_as_fifth_stage_in_frozen_order(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "aggregated_p48.json"
            target.write_text(
                json.dumps(
                    {
                        "sizes": [185, 265],
                        "independent_of_training_source": True,
                        "P4_S_prime": [0.0042, 0.0031],
                        "covariance_P4_S_prime": [[1e-8, 0.0], [0.0, 1.2e-8]],
                    }
                ),
                encoding="utf-8",
            )
            result = score(
                primary_score(), self.x17_artifact, self.p48_artifact, target
            )
        stage = result["stages"][4]
        self.assertEqual(stage["order"], 5)
        self.assertEqual(stage["status"], "SCORED_FROZEN_NO_REFIT")
        self.assertEqual(
            [row["name"] for row in stage["score"]["results"]],
            [
                "pure_power_baseline",
                "rank2_log_primary_correction",
                "analytic_inverse_N_competitor",
                "zero_effect",
            ],
        )


if __name__ == "__main__":
    unittest.main()
