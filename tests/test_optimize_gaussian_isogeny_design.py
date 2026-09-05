#!/usr/bin/env python3


from __future__ import annotations
import importlib.util
import sys
import unittest
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SPEC = importlib.util.spec_from_file_location(
    "optimize_gaussian_isogeny_design",
    ROOT / "scripts" / "optimize_gaussian_isogeny_design.py",
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class OptimalGaussianIsogenyDesignTest(unittest.TestCase):
    def test_norm5_exact_harmonic_fingerprints(self) -> None:
        row = MODULE.selected_campaign()[2]
        self.assertEqual(
            Fraction(
                row["exact_harmonics"]["H4"]["angular_ratio"]["numerator"],
                row["exact_harmonics"]["H4"]["angular_ratio"]["denominator"],
            ),
            Fraction(-14, 25),
        )
        self.assertEqual(
            Fraction(
                row["exact_harmonics"]["H12"]["angular_ratio"]["numerator"],
                row["exact_harmonics"]["H12"]["angular_ratio"]["denominator"],
            ),
            Fraction(23506, 15625),
        )

    def test_norm5_children_head_balanced_enumeration(self) -> None:
        rows = MODULE.enumerate_candidates(65)
        leading = {(row["parent"]["N"], row["child"]["N"]) for row in rows[:2]}
        self.assertEqual(leading, {(65, 325), (85, 425)})

    def test_frozen_maximin_allocation(self) -> None:
        result = MODULE.optimize_allocations(MODULE.selected_campaign())
        self.assertEqual(
            [row["replicas"] for row in result["allocations"]],
            [2_000_000_000, 600_000_000, 500_000_000, 500_000_000],
        )
        self.assertGreater(result["maximin_expected_chi_square"], 16.0)
        self.assertLess(result["used_billion_site_replicas"], 750.0)

    def test_norm4_is_explicit_noncyclic_benchmark(self) -> None:
        rows = MODULE.idealized_norm4_benchmark()
        self.assertEqual([row["child"]["N"] for row in rows], [260, 340])
        self.assertEqual(
            [row["child"]["smith_invariants"] for row in rows],
            [[2, 130], [2, 170]],
        )
        self.assertTrue(all(not row["child"]["cyclic"] for row in rows))
        for row in rows:
            self.assertEqual(
                Fraction(
                    row["exact_harmonics"]["H4"]["angular_ratio"]["numerator"],
                    row["exact_harmonics"]["H4"]["angular_ratio"]["denominator"],
                ),
                Fraction(1),
            )


if __name__ == "__main__":
    unittest.main()
