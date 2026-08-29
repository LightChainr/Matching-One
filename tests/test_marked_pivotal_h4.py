from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_marked_pivotal_h4 import read_batches, state  # noqa: E402


RESULT_ROOT = ROOT / "results/local-20260829/P100-marked-pivotal-h4"


class MarkedPivotalH4Tests(unittest.TestCase):
    def test_committed_exact_oracle_closes_russo_and_symmetries(self) -> None:
        payload = json.loads(
            (RESULT_ROOT / "analysis/exact-axis-l4-r1.json").read_text()
        )
        self.assertEqual(payload["fixed_root_configurations"], 2**15)
        self.assertEqual(float(payload["russo_control"]["primal_difference"]), 0.0)
        self.assertEqual(float(payload["russo_control"]["matching_difference"]), 0.0)
        self.assertEqual(
            payload["symmetry_violations"],
            {"rotation_C4": 0, "reflection": 0, "registry_pi_over_4": 0},
        )

    def test_committed_pilot_has_integer_identities_and_primary_marks(self) -> None:
        batches = read_batches(RESULT_ROOT / "raw/n65_r3_200k.batches.csv")
        self.assertEqual(len(batches["first"]), 100)
        for orientation in ("first", "second"):
            point = state(batches[orientation])
            self.assertGreater(point["mu4"], 0.0)
            self.assertGreater(point["a4"], 0.0)
            self.assertGreater(point["landing_acceptance_given_pivotal"], 0.5)

    def test_cpp_counter_smoke(self) -> None:
        compiler = shutil.which("c++")
        if compiler is None:
            self.skipTest("C++ compiler unavailable")
        with tempfile.TemporaryDirectory() as directory:
            binary = Path(directory) / "counter"
            prefix = Path(directory) / "pilot"
            subprocess.run(
                [compiler, "-std=c++17", "-O1", str(ROOT / "src/marked_pivotal_h4_mc.cpp"),
                 "-o", str(binary)],
                check=True,
            )
            subprocess.run(
                [str(binary), "--samples", "1000", "--batches", "10",
                 "--threads", "1", "--output-prefix", str(prefix)],
                check=True,
                capture_output=True,
                text=True,
            )
            batches = read_batches(Path(str(prefix) + ".batches.csv"))
            self.assertEqual(sum(row["samples"] for row in batches["first"]), 1000)


if __name__ == "__main__":
    unittest.main()
