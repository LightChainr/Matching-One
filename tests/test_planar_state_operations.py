from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from planar_state_operations import (  # noqa: E402
    detach_blocks,
    detach_rgs,
    join_cyclic_adjacent_blocks,
    join_cyclic_adjacent_rgs,
    validate_contract,
)


class PlanarStateOperationsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(
            (ROOT / "analysis" / "planar_state_operations_contract.json").read_text(
                encoding="utf-8"
            )
        )

    def test_checked_in_contract_exhausts_width_eight(self) -> None:
        result = validate_contract(self.contract)
        self.assertEqual(result["states_checked"], 2055)
        self.assertEqual(result["cases_per_operation"], 15521)
        self.assertEqual(result["combined_operation_cases"], 31042)
        self.assertTrue(result["independent_implementations_agree"])
        self.assertTrue(result["all_outputs_canonical_noncrossing"])
        self.assertTrue(result["all_operations_idempotent"])
        self.assertFalse(result["contains_transfer_matrix_result"])

    def test_detach_splits_a_point_and_is_idempotent(self) -> None:
        state = (0, 0, 1, 1, 0)
        expected = (0, 1, 2, 2, 1)
        self.assertEqual(detach_rgs(state, 0), expected)
        self.assertEqual(detach_blocks(state, 0), expected)
        self.assertEqual(detach_rgs(expected, 0), expected)

    def test_linear_and_wraparound_adjacent_joins(self) -> None:
        self.assertEqual(join_cyclic_adjacent_rgs((0, 1, 2, 1), 0), (0, 0, 1, 0))
        self.assertEqual(join_cyclic_adjacent_blocks((0, 1, 2, 1), 0), (0, 0, 1, 0))
        self.assertEqual(join_cyclic_adjacent_rgs((0, 1, 1, 2), 3), (0, 1, 1, 0))
        self.assertEqual(join_cyclic_adjacent_blocks((0, 1, 1, 2), 3), (0, 1, 1, 0))

    def test_operations_reject_crossing_states_and_bad_points(self) -> None:
        with self.assertRaisesRegex(ValueError, "crossing"):
            detach_rgs((0, 1, 0, 1), 0)
        for point in (-1, 4, True):
            with self.subTest(point=point):
                with self.assertRaises(ValueError):
                    join_cyclic_adjacent_rgs((0, 1, 1, 0), point)

    def test_contract_coverage_drift_fails_closed(self) -> None:
        changed = deepcopy(self.contract)
        changed["expected_cases_per_operation"] -= 1
        with self.assertRaisesRegex(ValueError, "operation coverage drifted"):
            validate_contract(changed)


if __name__ == "__main__":
    unittest.main()
