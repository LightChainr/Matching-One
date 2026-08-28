from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from c4_self_matching_exact import CHANNELS, enumerate_exact  # noqa: E402


class ExactC4SelfMatchingPolynomialCppTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        compiler = shutil.which("c++") or shutil.which("g++")
        if compiler is None:
            raise unittest.SkipTest("no C++17 compiler found")
        cls.tempdir = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.tempdir.name) / "exact_c4"
        subprocess.run(
            [
                compiler,
                "-O3",
                "-std=c++17",
                "-pthread",
                str(ROOT / "scripts" / "exact_c4_self_matching_polynomial.cpp"),
                "-o",
                str(cls.binary),
            ],
            check=True,
        )

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tempdir.cleanup()

    def run_cpp(self, threads: int) -> dict[str, object]:
        completed = subprocess.run(
            [
                str(self.binary),
                "--a",
                "3",
                "--b",
                "1",
                "--threads",
                str(threads),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        return json.loads(completed.stdout)

    def test_cpp_matches_python_all_n10_channels(self) -> None:
        cpp = self.run_cpp(2)
        python = enumerate_exact(3, 1)
        self.assertEqual(cpp["geometry"]["N"], 10)
        self.assertEqual(cpp["configurations"], 1024)
        for channel in CHANNELS:
            self.assertEqual(
                cpp["channels"][channel]["R_bernstein_integer_coefficients"],
                python["channels"][channel]["R_bernstein_integer_coefficients"],
            )
            self.assertEqual(
                cpp["channels"][channel]["M_bernstein_integer_coefficients"],
                python["channels"][channel]["M_bernstein_integer_coefficients"],
            )

    def test_thread_count_does_not_change_coefficients(self) -> None:
        single = self.run_cpp(1)
        multi = self.run_cpp(3)
        for channel in CHANNELS:
            self.assertEqual(single["channels"][channel], multi["channels"][channel])


if __name__ == "__main__":
    unittest.main()
