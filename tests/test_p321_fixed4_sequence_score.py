#!/usr/bin/env python3
"""Regression tests for the parameter-free P321 fixed-power diagnostic."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "p321_fixed4_sequence_score.py"
DATA = ROOT / "data" / "jacobsen_2015_square_site_cylinder.csv"


class P321FixedPowerSequenceScoreTests(unittest.TestCase):
    def test_committed_sequence_has_three_dyadic_ratios(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "score.json"
            subprocess.run(
                [sys.executable, str(SCRIPT), str(DATA), "--output", str(output)],
                cwd=ROOT,
                check=True,
                text=True,
                capture_output=True,
            )
            score = json.loads(output.read_text())
            rows = score["dyadic_diagnostics"]
            self.assertEqual([row["n"] for row in rows], [1, 2, 3, 4, 5])
            self.assertEqual(score["primary_hypothesis"]["dyadic_target"], "0.0625")
            self.assertFalse(score["primary_hypothesis"]["uses_fitted_p_infinity"])

    def test_exact_n_minus_four_sequence_hits_target(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            data = Path(directory) / "synthetic.csv"
            data.write_text(
                "n,value\n" + "\n".join(
                    f"{n},{1 + 3 / n**4:.17g}" for n in (1, 2, 4, 8)
                ) + "\n"
            )
            output = Path(directory) / "score.json"
            subprocess.run(
                [sys.executable, str(SCRIPT), str(data), "--output", str(output)],
                cwd=ROOT,
                check=True,
                text=True,
                capture_output=True,
            )
            score = json.loads(output.read_text())
            for row in score["dyadic_diagnostics"]:
                self.assertAlmostEqual(float(row["ratio"]), 1.0 / 16.0, places=13)


if __name__ == "__main__":
    unittest.main()
