#!/usr/bin/env python3
"""Compile the frontier C++ kernel and regress it against the Python oracle."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from exact_matching_polynomial import bernstein_counts, bernstein_to_power  # noqa: E402
from matched_torus_reference import axis_geometry  # noqa: E402


class ExactAxisMatchingPolynomialCppTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        compiler = shutil.which("c++") or shutil.which("clang++") or shutil.which("g++")
        if compiler is None:
            raise unittest.SkipTest("no C++ compiler")
        cls.tempdir = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.tempdir.name) / "exact-axis"
        subprocess.run(
            [compiler, "-std=c++17", "-O3", "-pthread",
             str(SCRIPTS / "exact_axis_matching_polynomial.cpp"), "-o", str(cls.binary)],
            check=True,
        )

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tempdir.cleanup()

    def test_l1_through_l4_match_python_exact_coefficients(self) -> None:
        for length in range(1, 5):
            completed = subprocess.run(
                [str(self.binary), "--L", str(length), "--threads", "2"],
                check=True,
                text=True,
                stdout=subprocess.PIPE,
            )
            result = json.loads(completed.stdout)
            expected_bernstein = bernstein_counts(axis_geometry(length))
            self.assertEqual(result["bernstein_counts"], expected_bernstein)
            self.assertEqual(
                result["power_coefficients_ascending"],
                bernstein_to_power(expected_bernstein),
            )


if __name__ == "__main__":
    unittest.main()
