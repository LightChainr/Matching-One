from __future__ import annotations

import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_p48_retrospective import read_histograms  # noqa: E402
from score_p200_n580_phaseA import (  # noqa: E402
    STATE_ORDER,
    intrinsic_statistics,
    load_prediction,
    render,
    score_models,
    state_from_statistics,
)
from score_p50_fullcurve_n290 import grouped  # noqa: E402


class P200N580ScorerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.prediction_path = (
            ROOT / "predictions" / "p200_n580_q2_jordan_score_input_20260829.json"
        )

    def test_existing_p50_fixture_reproduces_frozen_N145_state(self) -> None:
        path = ROOT / "results/server-20260829/P50-n145-n290-fullcurve/raw/n145_100m.hist.csv"
        records = read_histograms(path)
        point = state_from_statistics(intrinsic_statistics(grouped(records, 145)), 145)
        expected = [-0.008480909487795181, -0.011326674713668418,
                    0.4467428179653149, 1.4458006515540685]
        for actual, frozen in zip(point, expected):
            self.assertAlmostEqual(actual, frozen, places=13)

    def test_model_order_and_covariance_sum_are_frozen(self) -> None:
        prediction = load_prediction(self.prediction_path)
        point = [prediction["models"]["ordinary_q2"]["N580_state_prediction"][key]
                 for key in STATE_ORDER]
        target_covariance = [[0.01 if i == j else 0.0 for j in range(4)] for i in range(4)]
        rows = score_models(point, target_covariance, prediction)
        self.assertEqual([row["model"] for row in rows], ["ordinary_q2", "rank2_Jordan"])
        self.assertAlmostEqual(rows[0]["joint_GLS"]["chi_square"], 0.0)
        source = prediction["models"]["ordinary_q2"]["N580_state_prediction_covariance"]
        self.assertAlmostEqual(rows[0]["target_plus_prediction_covariance"][2][2],
                               0.01 + source[2][2])

    def test_synthetic_schema_fixture_runs_end_to_end(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            hist = root / "synthetic.hist.csv"
            moments = root / "synthetic.moments.csv"
            metadata = root / "synthetic.metadata.json"
            with hist.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, lineterminator="\n")
                writer.writerow(("n", "a", "b", "orientation", "batch", "samples", "kind", "k", "count"))
                for batch in range(4):
                    for orientation, rep, kminus, kplus in (
                        ("first", (24, 2), 276 + batch, 301 + batch),
                        ("second", (18, 16), 280 + 2 * batch, 297 + 2 * batch),
                    ):
                        writer.writerow((580, *rep, orientation, batch, 10, "minus", kminus, 10))
                        writer.writerow((580, *rep, orientation, batch, 10, "plus", kplus, 10))
            with moments.open("w", newline="", encoding="utf-8") as handle:
                fields = ("n", "a", "b", "orientation", "batch", "samples", "sum_kminus",
                          "sum_kplus", "sum_kminus2", "sum_kplus2", "sum_product",
                          "sum_gap", "sum_gap2")
                writer = csv.writer(handle, lineterminator="\n")
                writer.writerow(fields)
                for batch in range(4):
                    for orientation, rep, km, kp in (
                        ("first", (24, 2), 276 + batch, 301 + batch),
                        ("second", (18, 16), 280 + 2 * batch, 297 + 2 * batch),
                    ):
                        gap = kp - km
                        writer.writerow((580, *rep, orientation, batch, 10, 10 * km, 10 * kp,
                                         10 * km * km, 10 * kp * kp, 10 * km * kp,
                                         10 * gap, 10 * gap * gap))
            metadata.write_text(json.dumps({
                "git_commit": "synthetic-fixture", "seed": 7,
                "replica_counter_first": 1000, "replica_counter_last_exclusive": 1040,
                "samples_per_pair": 40, "batches": 4,
                "designs": [{
                    "N": 580, "first": [24, 2], "second": [18, 16],
                    "first_period_matrix": [[24, -2], [2, 24]],
                    "second_period_matrix": [[18, -16], [16, 18]],
                    "first_smith_invariants": [2, 290],
                    "second_smith_invariants": [2, 290],
                }],
            }), encoding="utf-8")
            payload = render(hist, moments, metadata, self.prediction_path)
            self.assertEqual(payload["status"], "scorer_frozen_before_N580_target_completion")
            self.assertEqual(payload["target"]["delete_one_batches"], 4)
            self.assertEqual([row["model"] for row in payload["models_in_frozen_order"]],
                             ["ordinary_q2", "rank2_Jordan"])


if __name__ == "__main__":
    unittest.main()
