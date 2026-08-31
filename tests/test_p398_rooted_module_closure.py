import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from p398_rooted_module_closure import build_result


class TestRootedModuleClosure(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result()

    def test_full_intersection_changes_only_after_coupling(self):
        first, second = self.result["candidates"]
        self.assertEqual((first["extended_dimension"], first["radical_dimension"], first["radical_Gram_skew_rank"]), (21, 6, 2))
        self.assertEqual(first["full_inherited_intersection"], "empty")
        self.assertEqual((second["extended_dimension"], second["radical_dimension"], second["radical_Gram_skew_rank"]), (23, 4, 0))
        self.assertEqual(second["full_inherited_intersection"], "unique")
        self.assertEqual(second["radical_leading_Gram_rank"], 4)

    def test_block_uniqueness_is_exact(self):
        proof = self.result["ordinary_block_uniqueness"]
        self.assertEqual(proof["ordinary_shifted_hom_rank"], 195)
        self.assertEqual(proof["ordinary_shifted_hom_dimension"], 1)
        self.assertEqual(proof["common_left_invariant_rank"], 13)
        self.assertEqual(proof["common_left_invariant_dimension"], 1)
        self.assertEqual(proof["ordinary_velocity_covariance_ranks"], [0] * 4)

    def test_residual_irrep_and_minimum_repair(self):
        first, second = self.result["candidates"]
        self.assertEqual(first["radical_C4_dimensions"], {"trivial": 3, "charge1_rational": 2, "charge2": 1})
        self.assertEqual(second["radical_C4_dimensions"], {"trivial": 3, "charge1_rational": 0, "charge2": 1})
        self.assertEqual(first["first_empty_certificate"]["affine_modulus_coefficient_rank"], 0)
        self.assertEqual(first["first_empty_certificate"]["augmented_rank"], 1)
        self.assertEqual(first["charge1_obstruction"]["leading_Gram_rank"], 2)
        self.assertEqual(second["charge1_obstruction"]["dimension"], 0)
        self.assertEqual(self.result["minimal_coupled_mark"]["additional_rational_dimension"], 2)
        self.assertFalse(self.result["minimal_coupled_mark"]["new_terminal_character"])

    def test_all_non_gram_gates_and_stronger_boundary(self):
        for row in self.result["candidates"]:
            for key, value in row["canonical_checks"].items():
                self.assertEqual(value, [0, 0, 0, 0] if isinstance(value, list) else 0, key)
            self.assertGreater(row["unprojected_G0_self_adjoint_residual_rank"], 0)
            self.assertGreater(row["unprojected_G1_self_adjoint_residual_rank"], 0)

    def test_active_reference_semantics_are_not_conflated(self):
        row = self.result["active_reference_distinction"]
        self.assertEqual((row["initial_dimension"], row["C4_closure_dimension"], row["join_closure_dimension"], row["join_detach_source_closure_dimension"]), (7, 7, 13, 14))
        self.assertEqual(row["minimal_witness"]["left_on_all_seven_inputs"], [0] * 7)
        self.assertEqual(row["minimal_witness"]["left_on_image"], 1)

    def test_committed_certificate_recomputes(self):
        stored = json.loads((ROOT / "results/p398-rooted-module-closure/latest.json").read_text())
        self.assertEqual(self.result, stored)


if __name__ == "__main__":
    unittest.main()
