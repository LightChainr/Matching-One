from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from z5_charged_ope_fusion import (  # noqa: E402
    build_artifact,
    conjugate_triple,
    determinant,
    neutral,
    primitive_neutral_triples,
)


class Z5ChargedOPEFusionTests(unittest.TestCase):
    def test_primitive_cubic_channels_and_conjugates(self) -> None:
        triples = primitive_neutral_triples()
        self.assertEqual(
            triples, [(1, 1, 3), (1, 2, 2), (2, 4, 4), (3, 3, 4)]
        )
        self.assertEqual(conjugate_triple((1, 1, 3)), (2, 4, 4))
        self.assertEqual(conjugate_triple((1, 2, 2)), (3, 3, 4))
        self.assertTrue(all(neutral(row) for row in triples))

    def test_cubic_phase_gauge_matrix_has_no_monomial_null_direction(self) -> None:
        matrix = [
            [2, 1, 0, 0],
            [0, 2, 1, 0],
            [1, 0, 0, 2],
            [0, 0, 2, 1],
        ]
        self.assertEqual(determinant(matrix), -15)

    def test_frozen_cubic_hecke_targets_are_distinct_and_on_unit_circle(self) -> None:
        artifact = build_artifact()
        rows = artifact["frozen_next_score"]["candidate_cubic_phases"]
        exact = set()
        for row in rows.values():
            self.assertIs(row["unit_modulus_check"], True)
            exact.add(tuple(row["exact_handed_ratio"].values()))
        self.assertEqual(len(exact), 3)

    def test_claim_boundary_keeps_exact_and_conditional_layers_separate(self) -> None:
        artifact = build_artifact()
        self.assertIn("determinant phase no-go", artifact["claim_boundary"]["exact"])
        self.assertIn("q_(3s)", artifact["claim_boundary"]["conditional_prediction"])
        self.assertEqual(
            artifact["frozen_next_score"]["joint_gls_degrees_of_freedom"], 4
        )


if __name__ == "__main__":
    unittest.main()
