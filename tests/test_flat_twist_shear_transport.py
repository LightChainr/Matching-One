from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_projective_birth_smoke import BirthCell  # noqa: E402
from score_flat_twist_shear_transport import (  # noqa: E402
    T_INVERSE_INTEGER,
    canonical_primitive,
    matmul2,
    score,
    transport_cell,
)


class FlatTwistShearTransportTests(unittest.TestCase):
    def test_period_and_primitive_line_transport(self) -> None:
        identity = [[8, -1], [1, 8]]
        self.assertEqual(matmul2(identity, T_INVERSE_INTEGER), [[8, -9], [1, 7]])
        self.assertEqual(canonical_primitive(-6, -4), (3, 2))
        cell = BirthCell(
            orientation="first", batch=0, samples=1,
            tau1=20, tau2=40, kind="LINE", ell_x=-3, ell_y=2, count=1,
        )
        self.assertEqual(transport_cell(cell), (20, 40, "LINE", 1, -2))

    def test_fresh_n65_smoke_closes_full_triplet(self) -> None:
        root = ROOT / "results/local-20260830/P337-flat-twist-shear-smoke"
        result = score(
            root / "shape_a.births.csv", root / "shape_a.metadata.json",
            root / "shape_b.births.csv", root / "shape_b.metadata.json",
            0.592746050790,
        )
        self.assertTrue(result["exact_gates"]["passed"])
        self.assertLess(result["exact_gates"]["max_batch_transport_residual"], 1e-15)
        joint = result["joint_identity_shear_estimate"]
        self.assertEqual(len(joint["order"]), 6)
        self.assertEqual(len(joint["covariance"]), 6)
        self.assertTrue(all(len(row) == 6 for row in joint["covariance"]))
        residual = result["transport_residual"]
        self.assertEqual(len(residual["mean"]), 3)
        self.assertLess(residual["max_absolute_mean"], 1e-16)
        self.assertIsNone(residual["joint_score"])


if __name__ == "__main__":
    unittest.main()
