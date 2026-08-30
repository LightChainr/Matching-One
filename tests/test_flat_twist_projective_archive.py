from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_flat_twist_projective_archive import (  # noqa: E402
    gate_residuals,
    projectivize,
    score,
    twist_characters,
    twist_sectors,
)


class FlatTwistProjectiveArchiveTests(unittest.TestCase):
    def test_primitive_lines_never_reduce_to_zero(self) -> None:
        self.assertEqual(projectivize(5, -2, 2), (1, 0))
        self.assertEqual(projectivize(5, -2, 3), (1, 2))
        self.assertEqual(projectivize(0, -1, 3), (0, 1))
        with self.assertRaises(ValueError):
            projectivize(6, 3, 3)

    def test_aggregate_trace_and_balanced_f3_character(self) -> None:
        state = {
            "P0": 0.2, "P1": 0.5, "P2": 0.3,
            "line_probabilities": {
                2: {(1, 0): 0.1, (1, 1): 0.3, (0, 1): 0.1},
                3: {(1, 0): 0.2, (1, 1): 0.1, (1, 2): 0.05, (0, 1): 0.15},
            },
            "raw_A4": 0j,
        }
        values = {"P0": 0.2, "P1": 0.5, "P2": 0.3}
        for q in (2, 3):
            sectors = twist_sectors(state, q)
            values[f"F{q}_S"] = sum(sectors.values())
            self.assertAlmostEqual(
                values[f"F{q}_S"], q * q * 0.2 + q * 0.5 + 0.3
            )
        characters = twist_characters(state)
        self.assertAlmostEqual(
            characters["F3_H4_axis_diag"],
            0.5 * (0.2 + 0.15 - 0.1 - 0.05),
        )
        self.assertTrue(all(abs(value) < 1e-14
                            for value in gate_residuals(state, values).values()))

    def test_n65_archive_full_covariance_and_sharpness(self) -> None:
        result_root = ROOT / "results/local-20260830/P334-projective-birth-N65-smoke"
        result = score(
            result_root / "n65_20k.births.csv",
            result_root / "n65_20k.metadata.json",
            0.592746050790,
        )
        self.assertTrue(result["sufficiency"]["archive_is_sufficient"])
        self.assertTrue(result["exact_gates"]["passed"])
        self.assertEqual(result["sufficiency"]["new_fields_required"], [])
        order = result["joint_estimate"]["order"]
        covariance = result["joint_estimate"]["covariance"]
        self.assertEqual(len(covariance), len(order))
        self.assertTrue(all(len(row) == len(order) for row in covariance))
        self.assertEqual(
            result["sharpness"]["same_H4_sector"]["winner"]["name"],
            "F3_H4_axis_diag",
        )
        self.assertEqual(
            result["sharpness"]["additional_twist_sectors"]["winner"]["name"],
            "F3_diagonal_odd",
        )


if __name__ == "__main__":
    unittest.main()

