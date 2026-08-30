from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import finite_abelian_twist_tomography as twist  # noqa: E402


class FiniteAbelianTwistTomographyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.certificate = twist.build_certificate()

    def test_saturated_annihilator_count_is_group_structure_independent(self) -> None:
        audit = self.certificate["machine_certificates"]["finite_abelian_groups"]
        self.assertTrue(audit["all_pass"])
        self.assertEqual(audit["row_count"], 30)
        for row in audit["rows"]:
            n = row["group_order"]
            self.assertEqual(row["counts"], {"rank0": n * n, "rank1": n, "rank2": 1})

    def test_order_two_three_vandermonde_inversion(self) -> None:
        probabilities = (Fraction(3, 17), Fraction(5, 17), Fraction(9, 17))
        trace_two = twist.aggregate_trace(2, probabilities)
        trace_three = twist.aggregate_trace(3, probabilities)
        self.assertEqual(
            twist.reconstruct_from_order_two_three(trace_two, trace_three),
            probabilities,
        )

    def test_each_prime_projective_line_has_q_minus_one_twists(self) -> None:
        rows = self.certificate["machine_certificates"]["prime_projective_orbits"]
        for row in rows:
            self.assertEqual(row["projective_line_count"], row["prime"] + 1)
            self.assertTrue(row["every_line_has_q_minus_1_twists"])

    def test_nonzero_twists_recover_modular_line_weights(self) -> None:
        rows = self.certificate["machine_certificates"]["prime_line_tomography"]
        for row in rows:
            self.assertTrue(row["aggregate_pass"])
            self.assertTrue(row["line_recovery_pass"])

    def test_checked_in_artifacts_reproduce(self) -> None:
        expected_json = json.loads(
            (ROOT / "results/finite-abelian-twist-tomography/latest.json").read_text(
                encoding="utf-8"
            )
        )
        expected_markdown = (
            ROOT / "results/finite-abelian-twist-tomography/latest.md"
        ).read_text(encoding="utf-8")
        self.assertEqual(expected_json, self.certificate)
        self.assertEqual(expected_markdown, twist.render_markdown(self.certificate) + "\n")


if __name__ == "__main__":
    unittest.main()
