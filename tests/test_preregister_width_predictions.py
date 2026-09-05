#!/usr/bin/env python3
"""Regression tests for the frozen width-22--24 prediction artifact."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "preregister_width_predictions.py"
DATA_RELATIVE = Path("data/jacobsen_2015_square_site_cylinder.csv")
FROZEN = ROOT / "predictions" / "polynomial_widths_22_24.yaml"
QUARANTINE = ROOT / "analysis" / "rational_stage_b_quarantine_manifest.json"

# SHA-256 of the frozen artifact with its input_sha256 line removed, i.e. of every
# preregistered prediction, interval and model coefficient.  It is pinned separately
# from the whole-file digest so that a correction to the *input table* -- which must
# move input_sha256 -- cannot quietly move a prediction with it.  This value predates
# the 2026-09-04 transcription correction and must survive any future one.
PREDICTION_BODY_SHA256 = "8fa0dae0973794c10a900e45c734957f22b31ddd7ba43ed2032fea3bb2acfe69"


def prediction_body_sha256(text: str) -> str:
    kept = [line for line in text.splitlines(keepends=True) if not line.startswith("input_sha256:")]
    return hashlib.sha256("".join(kept).encode("utf-8")).hexdigest()


class PreregisteredPredictionTests(unittest.TestCase):
    def test_script_reproduces_frozen_artifact_byte_for_byte(self) -> None:
        """A maintenance change must not silently alter the preregistration."""
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "prediction.yaml"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    str(DATA_RELATIVE),
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

    def test_predictions_survive_input_table_corrections(self) -> None:
        """Correcting the training table must never move a preregistered number.

        The pre-reveal barrier for widths 22--24 rests on the predictions being
        fixed before the targets are known.  A transcription fix to the input CSV
        legitimately changes ``input_sha256``; it must change nothing else.
        """
        text = FROZEN.read_text(encoding="utf-8")
        self.assertEqual(prediction_body_sha256(text), PREDICTION_BODY_SHA256)

    def test_quarantine_manifest_pins_the_current_frozen_artifact(self) -> None:
        manifest = json.loads(QUARANTINE.read_text(encoding="utf-8"))
        digest = hashlib.sha256(FROZEN.read_bytes()).hexdigest()
        self.assertEqual(manifest["frozen_prediction"]["sha256"], digest)
        correction = manifest["transcription_correction"]
        self.assertNotEqual(
            correction["frozen_prediction_sha256_before"],
            manifest["frozen_prediction"]["sha256"],
        )
        self.assertNotEqual(
            correction["training_source_sha256_before"],
            manifest["training_source"]["sha256"],
        )


if __name__ == "__main__":
    unittest.main()
