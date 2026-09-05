from __future__ import annotations

import copy
import json
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import jordan_nonseparation as oracle  # noqa: E402


class JordanNonseparationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.certificate = oracle.build_certificate()

    def test_exact_certificate_verifies(self) -> None:
        result = oracle.verify_certificate(self.certificate)
        self.assertEqual(result["ordinary_witnesses"], 6)
        self.assertEqual(result["context_error_inequalities"], 378)
        self.assertEqual(result["floating_point_operations"], 0)

    def test_true_jordan_endpoint_is_not_the_identity(self) -> None:
        a = oracle.transfer(F(1), F(0))
        n = oracle.difference(a, oracle.IDENTITY)
        self.assertEqual(oracle.multiply(n, n), ((F(0), F(0)), (F(0), F(0))))
        self.assertEqual(oracle.norm_inf(n), 1)

    def test_common_diagonalizer_for_two_cocycles(self) -> None:
        witness = self.certificate["witnesses"][0]
        s = oracle.matrix(witness["common_diagonalizer"])
        inverse = oracle.matrix(witness["diagonalizer_inverse"])
        for key in ("a", "b"):
            a = oracle.matrix(witness["ordinary_generators"][key])
            diagonal = oracle.multiply(oracle.multiply(inverse, a), s)
            self.assertEqual((diagonal[0][1], diagonal[1][0]), (0, 0))
            self.assertNotEqual(diagonal[0][0], diagonal[1][1])

    def test_composition_and_commutation_hold_exactly(self) -> None:
        generators = {"a": oracle.transfer(F(1), F(1, 100)),
                      "b": oracle.transfer(F(3), F(1, 100))}
        for u in oracle.words(("a", "b"), 2):
            for v in oracle.words(("a", "b"), 2):
                uv = oracle.word_matrix(u + v, generators)
                self.assertEqual(uv, oracle.multiply(oracle.word_matrix(u, generators),
                                                     oracle.word_matrix(v, generators)))
                self.assertEqual(uv, oracle.word_matrix(v + u, generators))

    def test_jordan_word_has_additive_cocycle(self) -> None:
        generators = {"a": oracle.transfer(F(1), F(0)), "b": oracle.transfer(F(3), F(0))}
        self.assertEqual(oracle.word_matrix(("a", "b", "b"), generators),
                         ((F(1), F(7)), (F(0), F(1))))

    def test_empty_context_has_zero_error_bound(self) -> None:
        self.assertEqual(oracle.word_matrix((), {}), oracle.IDENTITY)
        self.assertEqual(oracle.word_error_bound((), {}), 0)

    def test_worst_word_has_independent_closed_form(self) -> None:
        for epsilon in (F(1, 10), F(1, 100), F(1, 1000)):
            a = {"b": oracle.transfer(F(3), epsilon)}
            target = {"b": oracle.transfer(F(3), F(0))}
            delta = epsilon ** 2
            error = oracle.norm_inf(oracle.difference(oracle.word_matrix(("b",) * 5, a),
                                                      oracle.word_matrix(("b",) * 5, target)))
            self.assertEqual(error, 360 * delta + 648 * delta ** 2)

    def test_approximation_improves_while_conditioning_worsens(self) -> None:
        rows = self.certificate["witnesses"]
        for earlier, later in zip(rows, rows[1:]):
            self.assertLess(F(later["max_context_error_inf"]), F(earlier["max_context_error_inf"]))
            self.assertGreater(F(later["diagonalizer_condition_inf"]), F(earlier["diagonalizer_condition_inf"]))
        self.assertLess(F(rows[-1]["max_context_error_inf"]), F(1, 10 ** 9))

    def test_corrupted_generator_is_rejected(self) -> None:
        bad = copy.deepcopy(self.certificate)
        bad["witnesses"][0]["ordinary_generators"]["a"][1][0] = "0"
        with self.assertRaises(ValueError):
            oracle.verify_certificate(bad)

    def test_corrupted_diagonalizer_is_rejected(self) -> None:
        bad = copy.deepcopy(self.certificate)
        bad["witnesses"][1]["common_diagonalizer"][1][0] = "1"
        with self.assertRaises(ValueError):
            oracle.verify_certificate(bad)

    def test_corrupted_error_claim_is_rejected(self) -> None:
        bad = copy.deepcopy(self.certificate)
        bad["witnesses"][0]["max_context_error_inf"] = "0"
        with self.assertRaises(ValueError):
            oracle.verify_certificate(bad)

    def test_missing_witness_is_rejected(self) -> None:
        bad = copy.deepcopy(self.certificate)
        bad["witnesses"].pop()
        with self.assertRaises(ValueError):
            oracle.verify_certificate(bad)

    def test_float_certificate_entries_are_rejected(self) -> None:
        with self.assertRaises(TypeError):
            oracle.matrix(((1, 1.0), (0, 1)))

    def test_bad_shape_and_depth_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            oracle.matrix(((1, 0, 1), (0, 1, 1)))
        with self.assertRaises(ValueError):
            oracle.words(("a",), -1)

    def test_checked_in_certificate_reproduces(self) -> None:
        stored = json.loads((ROOT / "results/jordan-nonseparation/latest.json").read_text())
        self.assertEqual(stored, self.certificate)


if __name__ == "__main__":
    unittest.main()
