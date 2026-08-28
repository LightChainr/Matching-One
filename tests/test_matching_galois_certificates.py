from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from certify_matching_galois_groups import CASES, certify_case  # noqa: E402


class MatchingGaloisCertificateTests(unittest.TestCase):
    def test_full_symmetric_groups(self) -> None:
        expected = {
            ("axis", 3): "S_9",
            ("axis", 4): "S_16",
            ("axis", 5): "S_25",
            ("diamond", 2): "S_8",
            ("diamond", 3): "S_18",
            ("gaussian-3-1", 0): "S_10",
        }
        observed = {}
        for case in CASES:
            row = certify_case(case)
            key = (row["geometry"], row["L"])
            observed[key] = row["galois_group"]
            self.assertTrue(row["irreducible_over_Q"])
            self.assertTrue(row["transitive_galois_action"])
            self.assertTrue(row["factorization_squarefree"])
            self.assertTrue(row["primitive"])
            self.assertTrue(row["contains_transposition"])
            self.assertGreater(
                row["large_prime_cycle"],
                row["maximum_nontrivial_block_size_or_count"],
            )
            self.assertTrue(
                all(factor["irreducible_mod_prime"] for factor in row["factors"])
            )
        self.assertEqual(observed, expected)

    def test_cycle_types_isolate_transposition(self) -> None:
        for case in CASES:
            row = certify_case(case)
            self.assertEqual(row["factor_degrees"].count(2), 1)
            self.assertEqual(row["transposition_power"] % 2, 1)
            self.assertIn(row["large_prime_cycle"], row["factor_degrees"])


if __name__ == "__main__":
    unittest.main()
