from __future__ import annotations

import importlib.util
import json
from fractions import Fraction
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "p262_confluent_potts_projectors.py"
SPEC = importlib.util.spec_from_file_location("p262_projectors", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class ConfluentPottsProjectorTests(unittest.TestCase):
    def test_partition_gram_factorization(self) -> None:
        for n in (2, 4):
            for q in (4, 5, 6):
                direct = MODULE.determinant(MODULE.connectivity_gram(n, q))
                self.assertEqual(direct, MODULE.expected_gram_determinant(n, q))

    def test_Q1_connectivity_rank_loss(self) -> None:
        self.assertEqual(MODULE.matrix_rank(MODULE.connectivity_gram(2, 1)), 1)
        self.assertEqual(MODULE.matrix_rank(MODULE.connectivity_gram(4, 1)), 1)

    def test_pair_projectors_are_exact_orthogonal_idempotents(self) -> None:
        for q in (4, 5, 6, 7):
            projectors = MODULE.unordered_pair_projectors(q)
            self.assertTrue(MODULE.projectors_are_orthogonal(projectors.values()))
            self.assertEqual(MODULE.trace(projectors["singlet"]), 1)
            self.assertEqual(MODULE.trace(projectors["standard"]), q - 1)
            self.assertEqual(MODULE.trace(projectors["two_row_2"]), Fraction(q * (q - 3), 2))

    def test_projector_sum_is_identity(self) -> None:
        projectors = MODULE.unordered_pair_projectors(6)
        total = MODULE.add(*projectors.values())
        self.assertEqual(total, MODULE.identity(15))

    def test_machine_artifact_is_reproducible(self) -> None:
        expected = json.loads((ROOT / "predictions" / "p262_confluent_potts_projectors_20260829.json").read_text())
        self.assertEqual(MODULE.analyze(), expected)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "artifact.json"
            path.write_text(json.dumps(MODULE.analyze(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(path.read_bytes(), (ROOT / "predictions" / "p262_confluent_potts_projectors_20260829.json").read_bytes())


if __name__ == "__main__":
    unittest.main()
