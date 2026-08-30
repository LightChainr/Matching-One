from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from square_bond_primitive_norm2 import (  # noqa: E402
    LINEAGES,
    PRODUCTION_DESIGNS,
    ROTATE_ONE_PLUS_I,
    build_result,
    main,
    matrix_multiply,
    run_production_batches,
    seed_blocks,
    validate_designs,
)


class SquareBondPrimitiveNorm2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.rows, cls.blocks = run_production_batches(
            samples_per_design=160, batches=4, seed=41, workers=1
        )
        cls.result = build_result(
            cls.rows,
            cls.blocks,
            samples_per_design=160,
            batches=4,
            seed=41,
            dps=50,
        )

    def test_exact_norm2_geometry_map(self) -> None:
        validate_designs()
        expected_orders = [(30, 60), (56, 112)]
        for lineage, expected in zip(LINEAGES, expected_orders):
            parent = lineage["parent"][1]
            child = lineage["child"][1]
            self.assertEqual(matrix_multiply(ROTATE_ONE_PLUS_I, parent), child)
            determinants = []
            for matrix in (parent, child):
                determinants.append(
                    abs(
                        matrix[0][0] * matrix[1][1]
                        - matrix[0][1] * matrix[1][0]
                    )
                )
            self.assertEqual(tuple(determinants), expected)

    def test_four_independent_seed_blocks_and_sufficient_statistics(self) -> None:
        blocks = seed_blocks(41)
        self.assertEqual(len(blocks), 4)
        self.assertEqual(len(set(blocks.values())), 4)
        self.assertEqual(len(self.rows), 4 * 4)
        for row in self.rows:
            self.assertEqual(sum(row.counts.values()), row.samples)
            self.assertEqual(row.counts["invariant_failure"], 0)

    def test_parallel_worker_count_does_not_change_streams(self) -> None:
        serial, serial_blocks = run_production_batches(
            samples_per_design=80, batches=4, seed=73, workers=1
        )
        parallel, parallel_blocks = run_production_batches(
            samples_per_design=80, batches=4, seed=73, workers=2
        )
        self.assertEqual(serial_blocks, parallel_blocks)
        self.assertEqual(serial, parallel)

    def test_full_covariance_and_null_algebra(self) -> None:
        self.assertFalse(
            self.result["sampling_contract"]["parent_child_common_random_numbers"]
        )
        self.assertEqual(len(self.result["designs"]), len(PRODUCTION_DESIGNS))
        for design in self.result["designs"]:
            covariance = design["contrast_covariance_of_mean"]
            self.assertEqual((len(covariance), len(covariance[0])), (3, 3))
        for lineage in self.result["lineages"]:
            self.assertEqual(len(lineage["six_coordinate_covariance_of_mean"]), 6)
            cp = lineage["C_parent"]
            cc = lineage["C_child"]
            self.assertAlmostEqual(
                lineage["null_scores"]["H4_null_2Cchild_plus_Cparent"]["value"],
                2 * cc + cp,
            )
            self.assertAlmostEqual(
                lineage["null_scores"]["H8_null_2Cchild_minus_Cparent"]["value"],
                2 * cc - cp,
            )

    def test_cli_writes_json_and_batch_counts(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / "tiny"
            self.assertEqual(
                main(
                    [
                        "--samples-per-design",
                        "80",
                        "--batches",
                        "4",
                        "--workers",
                        "1",
                        "--dps",
                        "40",
                        "--output-prefix",
                        str(prefix),
                    ]
                ),
                0,
            )
            payload = json.loads(Path(str(prefix) + ".json").read_text())
            self.assertEqual(payload["samples_per_design"], 80)
            batch_lines = Path(str(prefix) + ".batches.csv").read_text().splitlines()
            self.assertEqual(len(batch_lines), 1 + 4 * 4)


if __name__ == "__main__":
    unittest.main()
