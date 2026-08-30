import json
from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from projection_leakage_commutator import (  # noqa: E402
    add,
    analyze,
    build_artifact,
    commutator,
    frobenius_norm_squared,
    identity,
    matrix,
    multiply,
    product,
    square,
    subtract,
    validate_artifact,
    witness_matrices,
)


class ProjectionLeakageCommutatorTest(unittest.TestCase):
    def test_committed_certificate_reproduces_exactly(self):
        path = ROOT / "analysis" / "projection_leakage_commutator_certificate.json"
        artifact = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(artifact, build_artifact())
        self.assertEqual(validate_artifact(artifact)["dimension"], 3)

    def test_full_translations_commute_and_are_involutions(self):
        u, v, _ = witness_matrices()
        zero = subtract(identity(3), identity(3))
        self.assertEqual(commutator(u, v), zero)
        self.assertEqual(multiply(u, u), identity(3))
        self.assertEqual(multiply(v, v), identity(3))

    def test_projection_is_exactly_idempotent(self):
        _, _, p = witness_matrices()
        self.assertEqual(multiply(p, p), p)

    def test_compressed_commutator_matches_frozen_matrix(self):
        u, v, p = witness_matrices()
        report = analyze(u, v, p)
        self.assertEqual(
            report["compressed_commutator"],
            [["0", "4/9", "4/9"], ["-4/9", "0", "4/9"], ["-4/9", "-4/9", "0"]],
        )
        self.assertEqual(report["compressed_commutator_frobenius_norm_squared"], "32/27")

    def test_nonzero_compressed_commutator_is_all_leakage(self):
        witness = build_artifact()["witness"]
        self.assertEqual(witness["full_commutator"], [["0"] * 3 for _ in range(3)])
        self.assertEqual(witness["intrinsic_projected_commutator"], [["0"] * 3 for _ in range(3)])
        self.assertEqual(witness["compressed_commutator"], witness["leakage_PVQUP_minus_PUQVP"])

    def test_commuting_projection_zero_control(self):
        control = build_artifact()["commuting_projection_zero_control"]
        self.assertEqual(control["compressed_commutator_frobenius_norm_squared"], "0")
        self.assertEqual(control["leakage_PVQUP_minus_PUQVP"], [["0"] * 3 for _ in range(3)])

    def test_matrix_helpers_close_exactly(self):
        a = matrix(((1, 2), (3, 4)))
        b = matrix(((4, 3), (2, 1)))
        self.assertEqual(add(a, b), matrix(((5, 5), (5, 5))))
        self.assertEqual(subtract(a, b), matrix(((-3, -1), (1, 3))))
        self.assertEqual(product(a, identity(2)), a)
        self.assertEqual(frobenius_norm_squared(a), Fraction(30))

    def test_float_and_boolean_values_fail_closed(self):
        with self.assertRaises(TypeError):
            matrix(((1.0, 0), (0, 1)))
        with self.assertRaises(TypeError):
            matrix(((True, 0), (0, 1)))

    def test_shape_and_projection_errors_fail_closed(self):
        with self.assertRaises(ValueError):
            matrix(((1, 2), (3,)))
        with self.assertRaises(ValueError):
            square(((1, 2, 3), (4, 5, 6)))
        u, v, _ = witness_matrices()
        with self.assertRaises(ValueError):
            analyze(u, v, identity(2))
        with self.assertRaises(ValueError):
            analyze(u, v, ((1, 0, 0), (0, 1, 0), (0, 0, "1/2")))

    def test_tampered_certificate_fails_validation(self):
        artifact = build_artifact()
        artifact["witness"]["compressed_commutator_frobenius_norm_squared"] = "0"
        with self.assertRaises(ValueError):
            validate_artifact(artifact)


if __name__ == "__main__":
    unittest.main()
