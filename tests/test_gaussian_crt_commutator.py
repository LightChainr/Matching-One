from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gaussian_crt_commutator import render  # noqa: E402


class GaussianCrtCommutatorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = render()
        cls.frozen = json.loads(
            (
                ROOT
                / "results"
                / "exact-cover-character-oracles"
                / "n650_gaussian_crt_commutator.json"
            ).read_text(encoding="utf-8")
        )

    def test_relative_smith_and_crt_square(self) -> None:
        factors = self.payload["factors"]
        self.assertEqual(factors["alpha_norm2"]["smith_invariants"], [1, 2])
        self.assertEqual(factors["beta_norm5"]["smith_invariants"], [1, 5])
        self.assertEqual(factors["gamma_product"]["smith_invariants"], [1, 10])
        self.assertEqual(factors["bezout_certificate"]["value"], [1, 0])
        self.assertTrue(self.payload["quotient_square"]["cartesian_and_cocartesian"])

    def test_n650_paths_have_expected_periods_and_smith_classes(self) -> None:
        paths = self.payload["n650_lattice_paths"]
        first, second = paths
        self.assertEqual(first["source"]["gaussian_period"], "8+i")
        self.assertEqual(
            first["path_B_norm2_then_norm5"]["intermediate"]["gaussian_period"],
            "7+9i",
        )
        self.assertEqual(
            first["path_A_norm5_then_norm2"]["intermediate"]["gaussian_period"],
            "17-6i",
        )
        self.assertEqual(first["path_A_norm5_then_norm2"]["final"]["coordinates"], [23, 11])
        self.assertEqual(second["path_A_norm5_then_norm2"]["final"]["coordinates"], [17, 19])
        for path in paths:
            self.assertEqual(path["source"]["smith_invariants"], [1, 65])
            self.assertEqual(
                path["path_B_norm2_then_norm5"]["intermediate"]["smith_invariants"],
                [1, 130],
            )
            self.assertEqual(
                path["path_A_norm5_then_norm2"]["intermediate"]["smith_invariants"],
                [1, 325],
            )
            self.assertEqual(
                path["path_A_norm5_then_norm2"]["final"]["smith_invariants"],
                [1, 650],
            )

    def test_linear_and_character_commutators_are_exactly_zero(self) -> None:
        self.assertTrue(
            self.payload["linear_fiber_oracle"]
            ["E5_E2_equals_E2_E5_equals_full_average"]
        )
        self.assertEqual(
            self.payload["linear_fiber_oracle"]["ordered_linear_mark"]["prediction"],
            "identically_zero_for_every_f",
        )
        self.assertTrue(
            self.payload["character_projector_oracle"]["commutator_zero_exact"]
        )
        self.assertEqual(
            len(self.payload["character_projector_oracle"]["all_joint_targets"]),
            10,
        )

    def test_full_joins_commute_but_mixed_nonlinear_defect_survives(self) -> None:
        oracle = self.payload["join_oracle"]
        self.assertTrue(oracle["full_partition_join_commutator_zero"])
        self.assertEqual(
            oracle["rank_values_for_discrete_input"],
            {"h_base": 0, "h_join_K2": 5, "h_join_K5": 8, "h_join_K2_K5": 9},
        )
        self.assertEqual(oracle["mixed_rank_defect"], -4)

    def test_frozen_artifact_matches(self) -> None:
        self.assertEqual(self.payload, self.frozen)


if __name__ == "__main__":
    unittest.main()
