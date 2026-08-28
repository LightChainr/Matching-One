#!/usr/bin/env python3
"""Regression tests for the frozen width-22--24 prediction artifact."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "preregister_width_predictions.py"
DATA = ROOT / "data" / "jacobsen_2015_square_site_cylinder.csv"
FROZEN = ROOT / "predictions" / "polynomial_widths_22_24.yaml"


class PreregisteredPredictionTests(unittest.TestCase):
    def test_script_reproduces_frozen_artifact_byte_for_byte(self) -> None:
        """A maintenance change must not silently alter the preregistration."""
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "prediction.yaml"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    str(DATA),
                    "--output",
                    str(output),
                    "--dps",
                    "100",
                ],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=True,
            )
            self.assertIn(str(output), completed.stdout)
            self.assertEqual(output.read_bytes(), FROZEN.read_bytes())


if __name__ == "__main__":
    unittest.main()
