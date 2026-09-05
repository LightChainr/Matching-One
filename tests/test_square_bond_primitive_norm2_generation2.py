
from __future__ import annotations
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from square_bond_primitive_norm2_generation2 import (  # noqa: E402
    LINEAGES,
    build_result,
    run_batches,
    score_model,
    validate_designs,
)


class Generation2Tests(unittest.TestCase):
    def test_exact_child_matrices(self) -> None:
        validate_designs()
        self.assertEqual(LINEAGES[0]["child_matrix"], ((0, -10), (12, 6)))
        self.assertEqual(LINEAGES[1]["child_matrix"], ((0, -14), (16, 8)))

    def test_fixed_score(self) -> None:
        parent = {
            "contrasts": {"C_nontrivial_real": {"value": 2.0}},
            "contrast_covariance_of_mean": [[0.04]],
        }
        child = {
            "contrasts": {"C_nontrivial_real": {"value": -1.0}},
            "contrast_covariance_of_mean": [[0.01]],
        }
        score = score_model(parent, child, -0.5)
        self.assertEqual(score["residual_Cchild_minus_ratio_Cparent"], 0.0)
        self.assertEqual(score["z"], 0.0)

    def test_tiny_parallel_is_deterministic(self) -> None:
        serial, serial_blocks = run_batches(
            samples_per_design=200, batches=4, seed=91, workers=1
        )
        parallel, parallel_blocks = run_batches(
            samples_per_design=200, batches=4, seed=91, workers=2
        )
        self.assertEqual(serial_blocks, parallel_blocks)
        self.assertEqual(serial, parallel)

    def test_build_uses_declared_generation1_parents(self) -> None:
        rows, blocks = run_batches(
            samples_per_design=200, batches=4, seed=92, workers=1
        )
        parent_payload = json.loads(
            (ROOT / "results" / "server-20260829" / "P156-norm2-h4-h8" / "result.json")
            .read_text(encoding="utf-8")
        )
        result = build_result(
            rows, blocks, parent_payload,
            samples_per_design=200, batches=4, seed=92, dps=40,
        )
        self.assertEqual(len(result["lineages"]), 2)
        self.assertEqual(set(result["joint_fixed_model_scores"]), {
            "rank4_H4", "even_nonlocal_character", "quadratic_H4",
            "local_H8_bound_saturated",
        })


if __name__ == "__main__":
    unittest.main()
