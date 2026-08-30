from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from noncrossing_connectivity_codec import (  # noqa: E402
    canonical_rgs,
    catalan,
    independent_noncrossing_states,
    is_noncrossing_rgs,
    noncrossing_states,
    rank_state,
    unrank_state,
    validate_contract,
    validate_rgs,
)


class NoncrossingConnectivityCodecTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(
            (ROOT / "analysis" / "noncrossing_connectivity_codec_contract.json").read_text(
                encoding="utf-8"
            )
        )

    def test_checked_in_contract_cross_checks_two_enumerators(self) -> None:
        result = validate_contract(self.contract)
        self.assertEqual(result["maximum_checked_width"], 8)
        self.assertTrue(result["independent_enumerators_agree"])
        self.assertTrue(result["all_rank_round_trips_exact"])
        self.assertFalse(result["contains_transfer_matrix_result"])

    def test_counts_are_catalan_through_width_eight(self) -> None:
        for width in range(1, 9):
            with self.subTest(width=width):
                self.assertEqual(len(noncrossing_states(width)), catalan(width))
                self.assertEqual(noncrossing_states(width), independent_noncrossing_states(width))

    def test_rank_and_unrank_are_exact(self) -> None:
        for width in range(1, 8):
            for expected_rank, state in enumerate(noncrossing_states(width)):
                self.assertEqual(rank_state(state), expected_rank)
                self.assertEqual(unrank_state(width, expected_rank), state)

    def test_crossing_and_nested_examples_are_distinguished(self) -> None:
        self.assertFalse(is_noncrossing_rgs((0, 1, 0, 1)))
        self.assertTrue(is_noncrossing_rgs((0, 1, 1, 0)))
        self.assertTrue(is_noncrossing_rgs((0, 0, 1, 1)))
        with self.assertRaisesRegex(ValueError, "crossing"):
            rank_state((0, 1, 0, 1))

    def test_canonicalization_and_malformed_inputs(self) -> None:
        self.assertEqual(canonical_rgs((7, 7, -2, 9, -2)), (0, 0, 1, 2, 1))
        for state in ((1,), (0, 2), (0, -1), ()):
            with self.subTest(state=state):
                with self.assertRaises(ValueError):
                    validate_rgs(state)
        with self.assertRaisesRegex(ValueError, "out of range"):
            unrank_state(3, 5)

    def test_contract_drift_fails_closed(self) -> None:
        changed = deepcopy(self.contract)
        changed["expected_counts"]["6"] = 131
        with self.assertRaisesRegex(ValueError, "declared Catalan count drifted"):
            validate_contract(changed)


if __name__ == "__main__":
    unittest.main()
