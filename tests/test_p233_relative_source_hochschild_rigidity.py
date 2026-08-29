from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p233_relative_source_hochschild_rigidity import build_oracle  # noqa: E402


class RelativeSourceHochschildRigidityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = build_oracle()
        cls.complex = cls.result["exact_complex"]

    def test_exact_cochain_dimensions(self) -> None:
        self.assertEqual(
            self.complex["dimensions"], {"C1": 9, "C2": 27, "C3": 81}
        )

    def test_hochschild_differential_squares_to_zero(self) -> None:
        self.assertTrue(self.complex["delta2_delta1_is_zero"])

    def test_second_cohomology_vanishes(self) -> None:
        self.assertEqual(self.complex["rank_delta1"], 9)
        self.assertEqual(self.complex["rank_delta2"], 18)
        self.assertEqual(self.complex["dimension_Z2"], 9)
        self.assertEqual(self.complex["dimension_B2"], 9)
        self.assertEqual(self.complex["dimension_HH2"], 0)

    def test_explicit_contraction_kills_every_cocycle_basis_vector(self) -> None:
        self.assertEqual(self.complex["nullspace_basis_size"], 9)
        self.assertTrue(
            self.complex["every_Z2_basis_vector_contracts_to_a_coboundary"]
        )

    def test_claim_points_outside_the_scalar_three_sector_algebra(self) -> None:
        consequence = self.result["research_consequence"]
        self.assertIn("larger non-semisimple", consequence["issue_233"])
        self.assertIn("connectivity-resolving", consequence["next_exact_target"])

    def test_committed_artifact_is_reproducible(self) -> None:
        committed = json.loads(
            (
                ROOT
                / "results"
                / "exact-relative-source-hochschild-rigidity"
                / "latest.json"
            ).read_text()
        )
        self.assertEqual(committed, self.result)


if __name__ == "__main__":
    unittest.main()
