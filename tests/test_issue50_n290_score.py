import csv
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from score_issue50_n290 import (  # noqa: E402
    CHANNELS,
    TARGET_DELTA_M,
    read_batches,
    score,
    validate_metadata,
    write_outputs,
)


class Issue50ScoreTests(unittest.TestCase):
    def setUp(self):
        self.commit = "a" * 40
        self.metadata = {
            "engine": "same-N Gaussian site-orientation discovery",
            "git_commit": self.commit,
            "command": (
                "build/gaussian_orientation_mc --samples 400 --batches 4 "
                "--n 290 --p-ref 0.59274605079 --seed 23 --replica-offset 900 "
                "--git-commit {} --output-prefix out"
            ).format(self.commit),
            "samples_per_pair": 400,
            "batches": 4,
            "p_ref": 0.59274605079,
            "seed": 23,
            "replica_counter_first": 900,
            "replica_counter_last_exclusive": 1300,
            "channels": list(CHANNELS),
            "designs": [{"N": 290, "first": [13, 11], "second": [17, 1]}],
        }

    def write_batches(self, path):
        fields = [
            "n", "batch", "samples", "p_ref", "channel", "a1", "b1", "a2", "b2",
            "first_primal_sum", "first_matching_sum", "second_primal_sum", "second_matching_sum",
        ]
        # Either-channel matching-function differences are .10,.20,.00,.10.
        first_primal = [60, 70, 50, 60]
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            for batch in range(4):
                for channel in CHANNELS:
                    writer.writerow({
                        "n": 290, "batch": batch, "samples": 100,
                        "p_ref": 0.59274605079, "channel": channel,
                        "a1": 13, "b1": 11, "a2": 17, "b2": 1,
                        "first_primal_sum": first_primal[batch] if channel == "either" else 50,
                        "first_matching_sum": 40 if channel == "either" else 50,
                        "second_primal_sum": 55 if channel == "either" else 50,
                        "second_matching_sum": 45 if channel == "either" else 50,
                    })

    def test_fixed_score_and_zero_benchmark(self):
        run = validate_metadata(self.metadata)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "run.batches.csv"
            self.write_batches(path)
            result = score(read_batches(path, run), run)
            json_path = Path(directory) / "score.json"
            csv_path = Path(directory) / "score.csv"
            write_outputs(result, json_path, csv_path)
            self.assertEqual(json.loads(json_path.read_text())["N"], 290)
            with csv_path.open(newline="", encoding="utf-8") as handle:
                self.assertEqual(next(csv.DictReader(handle))["lineage_first"], "13+11i")
        self.assertAlmostEqual(result["child_delta_M"], 0.1)
        self.assertAlmostEqual(result["child_sampling_se"], math.sqrt(0.02 / 12.0))
        self.assertAlmostEqual(
            result["target_residual_child_minus_frozen"], 0.1 - TARGET_DELTA_M
        )
        self.assertAlmostEqual(
            result["zero_z"], result["child_delta_M"] / result["child_sampling_se"]
        )

    def test_rejects_short_commit_and_counter_mismatch(self):
        for field, value, message in (
            ("git_commit", "abc1234", "complete 40-hex"),
            ("replica_counter_last_exclusive", 1299, "counter range"),
        ):
            with self.subTest(field=field):
                metadata = dict(self.metadata)
                metadata[field] = value
                with self.assertRaisesRegex(ValueError, message):
                    validate_metadata(metadata)

    def test_rejects_reversed_lineage(self):
        metadata = dict(self.metadata)
        metadata["designs"] = [{"N": 290, "first": [17, 1], "second": [13, 11]}]
        with self.assertRaisesRegex(ValueError, "lineage order"):
            validate_metadata(metadata)


if __name__ == "__main__":
    unittest.main()
