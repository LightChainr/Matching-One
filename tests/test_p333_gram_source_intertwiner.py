from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from p333_gram_source_intertwiner import (  # noqa: E402
    build_result,
    first_jet_radical_gram,
    join_block_count,
    rref_solve,
)
from noncrossing_connectivity_codec import noncrossing_states  # noqa: E402


RESULT = ROOT / "results/p333-gram-source-intertwiner/latest.json"


class TestP333GramSourceIntertwiner(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = json.loads(RESULT.read_text(encoding="utf-8"))

    def test_join_block_count(self):
        self.assertEqual(join_block_count((0, 1, 2), (0, 1, 2)), 3)
        self.assertEqual(join_block_count((0, 0, 1), (0, 1, 1)), 1)
        self.assertEqual(join_block_count((0, 0, 1), (0, 0, 1)), 2)

    def test_small_gram_matches_full_partition_first_jet(self):
        self.assertEqual(
            first_jet_radical_gram(noncrossing_states(2)),
            [[Fraction(1)]],
        )
        gram3 = first_jet_radical_gram(noncrossing_states(3))
        self.assertEqual(len(gram3), 4)
        self.assertEqual(len(gram3[0]), 4)

    def test_rref_affine_basis(self):
        solved = rref_solve([[1, 1]], [1], 2)
        self.assertTrue(solved["consistent"])
        self.assertEqual(solved["dimension"], 1)
        self.assertEqual(solved["particular"], [Fraction(1), Fraction(0)])
        self.assertEqual(solved["nullspace"], [[Fraction(-1), Fraction(1)]])

    def test_frozen_widths_and_basis_count(self):
        self.assertEqual([row["width"] for row in self.result["widths"]], [2, 3, 4])
        for row in self.result["widths"]:
            final = row["stages"]["source_normalized"]
            if final["consistent"]:
                self.assertEqual(
                    len(row["final_parameterization"]["primitive_integer_tangent_basis"]),
                    final["affine_tangent_dimension"],
                )

    def test_d0_affine_dimensions_reproduced(self):
        self.assertEqual(
            [row["stages"]["affine_sigma"]["affine_tangent_dimension"] for row in self.result["widths"]],
            [2, 3, 9],
        )

    def test_gram_and_source_gate(self):
        rows = self.result["widths"]
        self.assertEqual(
            [row["stages"]["gram_self_adjoint"]["affine_tangent_dimension"] for row in rows],
            [1, 1, 5],
        )
        self.assertEqual(
            [row["decision"] for row in rows],
            ["canonical_unique", "empty_intersection", "empty_intersection"],
        )
        self.assertTrue(rows[0]["canonical_translation"]["is_identity"])
        self.assertFalse(rows[1]["canonical_translation"]["is_identity"])
        self.assertFalse(rows[2]["canonical_translation"]["is_identity"])

    def test_exact_source_inconsistency_witnesses(self):
        for row in self.result["widths"][1:]:
            witness = row["source_restriction_on_gram_moduli"]["inconsistency_witness"]
            self.assertEqual(witness["left_times_rhs"], 1)
            self.assertTrue(all(value == 0 for value in witness["left_times_parameter_matrix"]))

    def test_full_certificate_recomputes(self):
        self.assertEqual(build_result(), self.result)


if __name__ == "__main__":
    unittest.main()
