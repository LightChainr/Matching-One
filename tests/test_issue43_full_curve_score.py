import csv
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_issue43_full_curve import (  # noqa: E402
    DESIGNS,
    P_REF,
    analyze,
    read_histograms,
    reconstruct,
    tail,
    validate_metadata,
    validate_moments,
)


class Issue43FullCurveTests(unittest.TestCase):
    def test_threshold_cdf_reconstruction(self):
        p = 0.37
        self.assertAlmostEqual(tail([0, 0, 1, 0], 1, p), 3 * p * p - 2 * p ** 3)

    def write_run(self, root, n, commit="a" * 40):
        first, second = DESIGNS[n]
        prefix = root / "n{}".format(n)
        hist_path = Path(str(prefix) + ".hist.csv")
        moments_path = Path(str(prefix) + ".moments.csv")
        metadata_path = Path(str(prefix) + ".metadata.json")
        hist_fields = ["n", "a", "b", "orientation", "batch", "samples", "kind", "k", "count"]
        moment_fields = [
            "n", "a", "b", "orientation", "batch", "samples", "sum_kminus", "sum_kplus",
            "sum_kminus2", "sum_kplus2", "sum_product", "sum_gap", "sum_gap2",
        ]
        moment_rows = []
        with hist_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=hist_fields, lineterminator="\n")
            writer.writeheader()
            for batch in range(100):
                shift = batch % 7 - 3
                for orientation, rep in (("first", first), ("second", second)):
                    distributions = {
                        "minus": (45 - shift, 55 + shift) if orientation == "first" else (50, 50),
                        "plus": (50 + shift, 50 - shift) if orientation == "first" else (50, 50),
                    }
                    sums = {}
                    for kind, counts in distributions.items():
                        writer.writerow({"n": n, "a": rep[0], "b": rep[1], "orientation": orientation, "batch": batch, "samples": 100, "kind": kind, "k": 1, "count": counts[0]})
                        writer.writerow({"n": n, "a": rep[0], "b": rep[1], "orientation": orientation, "batch": batch, "samples": 100, "kind": kind, "k": n, "count": counts[1]})
                        sums[kind] = (counts[0] + n * counts[1], counts[0] + n * n * counts[1])
                    moment_rows.append({
                        "n": n, "a": rep[0], "b": rep[1], "orientation": orientation,
                        "batch": batch, "samples": 100,
                        "sum_kminus": sums["minus"][0], "sum_kplus": sums["plus"][0],
                        "sum_kminus2": sums["minus"][1], "sum_kplus2": sums["plus"][1],
                        "sum_product": 0, "sum_gap": sums["plus"][0] - sums["minus"][0],
                        "sum_gap2": 0,
                    })
        with moments_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=moment_fields, lineterminator="\n")
            writer.writeheader(); writer.writerows(moment_rows)
        metadata = {
            "engine": "same-N Gaussian threshold-rank Newman-Ziff", "git_commit": commit,
            "command": "engine --samples 10000 --batches 100 --n {} --seed 51 --replica-offset {} --git-commit {} --output-prefix out".format(n, n * 10000, commit),
            "samples_per_pair": 10000, "batches": 100, "seed": 51,
            "replica_counter_first": n * 10000,
            "replica_counter_last_exclusive": n * 10000 + 10000,
            "designs": [{"N": n, "first": list(first), "second": list(second)}],
        }
        metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
        return hist_path, moments_path, metadata_path

    def test_accepts_real_microcanonical_metadata_schema_without_p_ref(self):
        with tempfile.TemporaryDirectory() as directory:
            _, _, metadata_path = self.write_run(Path(directory), 185)
            metadata = json.loads(metadata_path.read_text())
            metadata.update({
                "generated_utc": "2026-08-28T11:19:38Z", "compiler": "10.3.1",
                "openmp": True, "threads_requested": 16,
                "rng": "counter-derived SplitMix64 stream plus unbiased Fisher-Yates",
                "coupling": "same cyclic permutation shared by same-N orientations",
                "channel": "rank-2 cross wrapping",
                "K_plus": "first black primal cross rank, 1-based",
                "K_minus": "first black rank after white matching cross is lost; N-r+1",
                "sparse_joint_histogram": False, "per_batch_joint_moments": True,
                "elapsed_seconds": 65.0,
                "histogram_csv": "n185.hist.csv", "moments_csv": "n185.moments.csv",
            })
            self.assertNotIn("p_ref", metadata)
            metadata_path.write_text(json.dumps(metadata))
            run = validate_metadata(metadata_path)
            self.assertEqual(run["scorer_p_ref"], P_REF)

    def test_reconstruction_and_correlated_source_score(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            runs = {}
            for n in (185, 265):
                hist, moments, metadata_path = self.write_run(root, n)
                metadata = validate_metadata(metadata_path)
                records = read_histograms(hist, metadata)
                validate_moments(moments, metadata, records)
                sectors = reconstruct(records)
                self.assertEqual(len(sectors["DeltaM"]), 100)
                self.assertTrue(all(math.isfinite(value) for value in sectors["DeltaS"]))
                runs[n] = {"metadata": metadata, "sectors": sectors}
            result = analyze(runs, ROOT / "predictions" / "two_spin4_heldout_20260828.yaml")
            covariance = result["scores"]["DeltaM"]["target_covariance"]
            self.assertGreater(covariance[0][1], 0.0)
            self.assertEqual(result["scores"]["DeltaM"]["target_df"], 2)

    def test_rejects_short_commit_and_wrong_design_order(self):
        with tempfile.TemporaryDirectory() as directory:
            _, _, metadata_path = self.write_run(Path(directory), 185, commit="abc1234")
            with self.assertRaisesRegex(ValueError, "full 40-hex"):
                validate_metadata(metadata_path)
            _, _, metadata_path = self.write_run(Path(directory), 265)
            metadata = json.loads(metadata_path.read_text())
            metadata["designs"][0]["first"], metadata["designs"][0]["second"] = metadata["designs"][0]["second"], metadata["designs"][0]["first"]
            metadata_path.write_text(json.dumps(metadata))
            with self.assertRaisesRegex(ValueError, "orientation order"):
                validate_metadata(metadata_path)


if __name__ == "__main__":
    unittest.main()
