from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from analyze_threshold_rank_orientation import read_histograms  # noqa: E402
from integer_period_torus import gaussian_integer_torus  # noqa: E402
from threshold_rank_nz import simulate  # noqa: E402


class ThresholdRankOrientationMCTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        compiler = os.environ.get("CXX") or shutil.which("g++") or shutil.which("clang++")
        if compiler is None:
            raise unittest.SkipTest("no C++17 compiler found")
        cls.temporary = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.temporary.name) / "threshold_rank_orientation_mc"
        command = [compiler, "-O2", "-std=c++17"]
        if sys.platform != "darwin" and "clang" not in Path(compiler).name:
            command.append("-fopenmp")
        command += [
            str(ROOT / "src" / "threshold_rank_orientation_mc.cpp"),
            "-o", str(cls.binary),
        ]
        subprocess.run(command, check=True, cwd=ROOT)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    def test_exact_self_test(self) -> None:
        completed = subprocess.run(
            [str(self.binary), "--self-test"],
            check=True, text=True, capture_output=True,
        )
        self.assertIn("N=5 all permutations", completed.stdout)
        self.assertIn("Python-compatible", completed.stdout)

    def test_frozen_issue43_designs_are_available(self) -> None:
        for n, first, second in (
            (185, [13, 4], [11, 8]),
            (265, [16, 3], [12, 11]),
        ):
            prefix = Path(self.temporary.name) / f"n{n}"
            subprocess.run(
                [
                    str(self.binary),
                    "--samples", "4", "--batches", "2", "--n", str(n),
                    "--seed", "43", "--threads", "1",
                    "--output-prefix", str(prefix),
                ],
                check=True,
            )
            metadata = json.loads(
                Path(str(prefix) + ".metadata.json").read_text(encoding="utf-8")
            )
            self.assertEqual(metadata["designs"][0]["N"], n)
            self.assertEqual(metadata["designs"][0]["first"], first)
            self.assertEqual(metadata["designs"][0]["second"], second)

    def test_frozen_norm5_designs_are_available(self) -> None:
        for n, first, second in (
            (325, [17, 6], [18, 1]),
            (425, [16, 13], [19, 8]),
        ):
            prefix = Path(self.temporary.name) / f"norm5_n{n}"
            subprocess.run(
                [
                    str(self.binary),
                    "--samples", "4", "--batches", "2", "--n", str(n),
                    "--seed", "57", "--threads", "1",
                    "--output-prefix", str(prefix),
                ],
                check=True,
            )
            metadata = json.loads(
                Path(str(prefix) + ".metadata.json").read_text(encoding="utf-8")
            )
            self.assertEqual(metadata["designs"][0]["N"], n)
            self.assertEqual(metadata["designs"][0]["first"], first)
            self.assertEqual(metadata["designs"][0]["second"], second)

    def test_cpp_matches_python_and_is_thread_reproducible(self) -> None:
        first_prefix = Path(self.temporary.name) / "first"
        second_prefix = Path(self.temporary.name) / "second"
        common = [
            "--samples", "80", "--batches", "4", "--n", "65",
            "--seed", "17", "--replica-offset", "5", "--threads", "1",
            "--git-commit", "test-sha",
        ]
        subprocess.run(
            [str(self.binary), *common, "--output-prefix", str(first_prefix)],
            check=True,
        )
        threaded = list(common)
        threaded[threaded.index("--threads") + 1] = "2"
        subprocess.run(
            [str(self.binary), *threaded, "--output-prefix", str(second_prefix)],
            check=True,
        )
        for suffix in (".hist.csv", ".moments.csv"):
            self.assertEqual(
                Path(str(first_prefix) + suffix).read_bytes(),
                Path(str(second_prefix) + suffix).read_bytes(),
            )

        records = read_histograms(Path(str(first_prefix) + ".hist.csv"))
        for orientation, a, b in (("first", 8, 1), ("second", 7, 4)):
            selected = [records[key] for key in sorted(records) if key[1] == orientation]
            cpp_minus = [0] * 66
            cpp_plus = [0] * 66
            for record in selected:
                for rank in range(66):
                    cpp_minus[rank] += record["minus"][rank]
                    cpp_plus[rank] += record["plus"][rank]
            python = simulate(
                gaussian_integer_torus(a, b), 80, seed=17, counter_start=5
            )
            self.assertEqual(cpp_minus, python.kminus)
            self.assertEqual(cpp_plus, python.kplus)

        metadata = json.loads(
            Path(str(first_prefix) + ".metadata.json").read_text(encoding="utf-8")
        )
        self.assertEqual(metadata["git_commit"], "test-sha")
        self.assertEqual(metadata["replica_counter_first"], 5)
        self.assertTrue(metadata["per_batch_joint_moments"])

        analysis_json = Path(self.temporary.name) / "analysis.json"
        analysis_csv = Path(self.temporary.name) / "analysis.csv"
        subprocess.run(
            [
                sys.executable,
                str(SCRIPTS / "analyze_threshold_rank_orientation.py"),
                "--histograms", str(first_prefix) + ".hist.csv",
                "--moments", str(first_prefix) + ".moments.csv",
                "--p", "0.591746050790", "--p", "0.592746050790",
                "--json", str(analysis_json), "--csv", str(analysis_csv),
            ],
            check=True,
        )
        payload = json.loads(analysis_json.read_text(encoding="utf-8"))["by_N"]["65"]
        self.assertEqual(payload["first_rep"], [8, 1])
        self.assertEqual(payload["second_rep"], [7, 4])
        self.assertEqual(len(payload["evaluations"]), 2)
        self.assertEqual(payload["amplitude_closure"]["delete_one_batches"], 4)
        self.assertIn("C_jackknife_se", payload["amplitude_closure"])
        with analysis_csv.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual({row["metric"] for row in rows}, {"M", "M_prime", "root"})


if __name__ == "__main__":
    unittest.main()