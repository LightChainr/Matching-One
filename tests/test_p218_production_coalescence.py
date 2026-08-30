from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_p218_production_coalescence import (  # noqa: E402
    fit_jordan,
    matrix_diagnostics,
)


class P218ProductionCoalescenceTests(unittest.TestCase):
    def test_exact_jordan_certificate(self) -> None:
        result = matrix_diagnostics([[1.0, 2.0], [0.0, 1.0]])
        self.assertEqual(result["relative_eigenvalue_gap"], 0.0)
        self.assertEqual(result["J2_over_J1"], 0.0)
        self.assertEqual(result["right_eigenvector_principal_angle_degrees"], 0.0)

    def test_jordan_fit_recovers_synthetic_map(self) -> None:
        parent = [[1.0, 0.0], [0.0, 1.0]]
        child = [[1.25, 0.75], [0.0, 1.25]]
        fitted, details = fit_jordan(parent, child)
        for i in range(2):
            for j in range(2):
                self.assertAlmostEqual(fitted[i][j], child[i][j], places=8)
        self.assertLess(details["source_fit_squared_residual"], 1e-15)

    def test_committed_artifact_has_required_joint_diagnostics(self) -> None:
        payload = json.loads(
            (ROOT / "results/p218-production-coalescence/latest.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(len(payload["plug_in_transfer_geometry"]), 3)
        for row in payload["plug_in_transfer_geometry"]:
            self.assertIn("relative_eigenvalue_gap", row)
            self.assertIn("right_eigenbasis_condition_2norm", row)
            self.assertIn("J2_over_J1", row)
        self.assertEqual(
            set(payload["heldout_model_table"]),
            {"normal_diagonalizable", "rank2_Jordan", "generic_2x2"},
        )


if __name__ == "__main__":
    unittest.main()
