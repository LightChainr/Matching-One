
from __future__ import annotations
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_issue43_full_curve_locked import (  # noqa: E402
    FROZEN_BATCHES,
    FROZEN_COUNTERS,
    FROZEN_SAMPLES,
    FROZEN_SEED,
    FROZEN_SOURCE_COMMIT,
    FROZEN_THREADS,
    validate_metadata,
    validate_moments,
)


def metadata(n: int) -> dict:
    designs = {
        185: ([13, 4], [11, 8]),
        265: ([16, 3], [12, 11]),
    }
    first, second = designs[n]
    counter_first, counter_last = FROZEN_COUNTERS[n]
    command = (
        "./threshold_rank_orientation_mc "
        f"--n {n} --samples {FROZEN_SAMPLES} --batches {FROZEN_BATCHES} "
        f"--seed {FROZEN_SEED} --replica-offset {counter_first} "
        f"--threads {FROZEN_THREADS} --git-commit {FROZEN_SOURCE_COMMIT} "
        "--output-prefix frozen"
    )
    return {
        "engine": "same-N Gaussian threshold-rank Newman-Ziff",
        "git_commit": FROZEN_SOURCE_COMMIT,
        "command": command,
        "compiler": "g++ test compiler",
        "openmp": True,
        "threads_requested": FROZEN_THREADS,
        "samples_per_pair": FROZEN_SAMPLES,
        "batches": FROZEN_BATCHES,
        "seed": FROZEN_SEED,
        "replica_counter_first": counter_first,
        "replica_counter_last_exclusive": counter_last,
        "designs": [{"N": n, "first": first, "second": second}],
    }


class Issue43ProductionLockTests(unittest.TestCase):
    def write_metadata(self, payload: dict) -> Path:
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        path = Path(directory.name) / "metadata.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        return path

    def test_exact_frozen_metadata_passes_for_both_targets(self) -> None:
        for n in (185, 265):
            with self.subTest(n=n):
                run = validate_metadata(self.write_metadata(metadata(n)))
                self.assertEqual(run["N"], n)
                self.assertEqual(run["commit"], FROZEN_SOURCE_COMMIT)

    def test_scientifically_relevant_freeze_fields_are_not_merely_self_consistent(self) -> None:
        mutations = [
            ("git_commit", "0" * 40),
            ("samples_per_pair", FROZEN_SAMPLES // 2),
            ("seed", FROZEN_SEED + 1),
            ("replica_counter_first", FROZEN_COUNTERS[185][0] + 1),
            ("replica_counter_last_exclusive", FROZEN_COUNTERS[185][1] + 1),
            ("threads_requested", FROZEN_THREADS - 1),
            ("openmp", False),
        ]
        for key, value in mutations:
            with self.subTest(key=key):
                payload = metadata(185)
                payload[key] = value
                # Keep the command/self-consistency fields synchronized when
                # appropriate: the strict lock should still reject the run.
                if key == "git_commit":
                    payload["command"] = payload["command"].replace(
                        FROZEN_SOURCE_COMMIT, str(value)
                    )
                elif key == "samples_per_pair":
                    payload["replica_counter_last_exclusive"] = (
                        payload["replica_counter_first"] + int(value)
                    )
                    payload["command"] = payload["command"].replace(
                        f"--samples {FROZEN_SAMPLES}", f"--samples {value}"
                    )
                elif key == "seed":
                    payload["command"] = payload["command"].replace(
                        f"--seed {FROZEN_SEED}", f"--seed {value}"
                    )
                elif key == "replica_counter_first":
                    payload["replica_counter_last_exclusive"] = int(value) + FROZEN_SAMPLES
                    payload["command"] = payload["command"].replace(
                        f"--replica-offset {FROZEN_COUNTERS[185][0]}",
                        f"--replica-offset {value}",
                    )
                with self.assertRaises(ValueError):
                    validate_metadata(self.write_metadata(payload))

    def test_command_thread_count_is_locked_separately_from_metadata(self) -> None:
        payload = metadata(185)
        payload["command"] = payload["command"].replace(
            f"--threads {FROZEN_THREADS}", f"--threads {FROZEN_THREADS - 1}"
        )
        with self.assertRaisesRegex(ValueError, "threads"):
            validate_metadata(self.write_metadata(payload))

    def test_joint_squared_gap_identity_is_checked(self) -> None:
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        path = Path(directory.name) / "moments.csv"
        fields = [
            "n", "a", "b", "orientation", "batch", "samples",
            "sum_kminus", "sum_kplus", "sum_kminus2", "sum_kplus2",
            "sum_product", "sum_gap", "sum_gap2",
        ]
        row = {
            "n": 5, "a": 2, "b": 1, "orientation": "first", "batch": 0,
            "samples": 10, "sum_kminus": 30, "sum_kplus": 40,
            "sum_kminus2": 90, "sum_kplus2": 160, "sum_product": 120,
            "sum_gap": 10, "sum_gap2": 10,
        }
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerow(row)

        minus = [0] * 6
        plus = [0] * 6
        minus[3] = 10
        plus[4] = 10
        records = {("first", 0): {"samples": 10, "minus": minus, "plus": plus}}
        run = {"N": 5, "first": (2, 1), "second": (2, 1)}
        validate_moments(path, run, records)

        row["sum_gap2"] = 11
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerow(row)
        with self.assertRaisesRegex(ValueError, "squared-gap"):
            validate_moments(path, run, records)


if __name__ == "__main__":
    unittest.main()
