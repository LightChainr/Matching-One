from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class P267EtopTauTopologyScoreTests(unittest.TestCase):
    def test_locked_score_and_primary_boundary(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "score.json"
            report = Path(directory) / "REPORT.md"
            subprocess.run([
                sys.executable, str(ROOT / "scripts" / "score_p267_etop_tau_topology_factorial.py"),
                "--raw-dir", str(ROOT / "results" / "p267-etop-tau-topology-factorial" / "raw"),
                "--json", str(output), "--report", str(report),
            ], check=True, cwd=ROOT)
            payload = json.loads(output.read_text())
            self.assertEqual(payload["field_order"], ["A_top", "E_top", "C", "W"])
            self.assertEqual(payload["primary_character_normalized_interaction"]["df"], 4)
            self.assertEqual(payload["promotion_gate"]["pilot_decision"], "extend_both_to_100k")
            self.assertEqual(payload["promotion_gate"]["final_samples_per_missing_cell"], 100000)
            self.assertIn("no field identity", payload["claim_boundary"])


if __name__ == "__main__":
    unittest.main()
