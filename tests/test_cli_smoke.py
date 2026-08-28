#!/usr/bin/env python3
"""Cheap smoke coverage for high-use research command-line tools."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]

# These are the high-use analysis/experiment entry points.  The smoke contract is
# intentionally small: the module must import, construct its CLI, and answer
# --help without starting a production computation.
SCRIPTS = (
    "scripts/finite_size_audit.py",
    "scripts/rational_finite_size_audit.py",
    "scripts/analyze_gaussian_orientation_mc.py",
    "scripts/analyze_gaussian_doubling.py",
    "scripts/analyze_threshold_ranks.py",
    "scripts/analyze_threshold_rank_orientation.py",
    "scripts/challenge_orientation_radial_models.py",
    "scripts/score_angular_root_amplitude.py",
    "scripts/score_issue50_n290.py",
    "scripts/exact_matching_polynomial.py",
    "scripts/matched_torus_reference.py",
    "scripts/preregister_width_predictions.py",
    "scripts/run_finite_size_grid.py",
    "scripts/summarize_finite_size_grid.py",
)


class CliSmokeTests(unittest.TestCase):
    def test_high_use_clis_answer_help(self) -> None:
        for relative in SCRIPTS:
            with self.subTest(script=relative):
                path = ROOT / relative
                self.assertTrue(path.is_file(), relative)
                result = subprocess.run(
                    [sys.executable, str(path), "--help"],
                    cwd=ROOT,
                    text=True,
                    capture_output=True,
                    timeout=15,
                    check=False,
                )
                self.assertEqual(
                    result.returncode,
                    0,
                    msg=f"{relative}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}",
                )
                self.assertTrue(
                    result.stdout.strip() or result.stderr.strip(),
                    msg=f"{relative} produced no help output",
                )


if __name__ == "__main__":
    unittest.main()
