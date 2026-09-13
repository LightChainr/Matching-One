"""Tripwire: wrapping-type collapse reproduces committed Bernstein integers."""
from __future__ import annotations

import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class TestWrappingTypeCensus(unittest.TestCase):
    def test_axis_l3_artifact_matches_committed_gate(self) -> None:
        path = ROOT / "results" / "wrapping-type-census" / "axis-L3.json"
        data = json.loads(path.read_text())
        self.assertEqual(data["tripwire"], "pass")
        self.assertEqual(
            data["collapsed"],
            [-1, -9, -36, -78, -90, -36, 36, 36, 9, 1],
        )
        self.assertEqual(data["tables"]["5"]["neither"]["both"], 45)
        self.assertEqual(data["tables"]["5"]["both"]["neither"], 9)
