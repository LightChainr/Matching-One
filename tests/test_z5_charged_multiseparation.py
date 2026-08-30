from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_p250_near_zero_reveal import audit  # noqa: E402
from score_z5_charged_multiseparation import read_batches, score  # noqa: E402
from z5_charged_multiseparation_mc import (  # noqa: E402
    FIELD_ORDER,
    SMOKE_CAP,
    exact_mapping_gate,
    hermitian_pair,
    run,
)


RESULT = ROOT / "results/local-20260830/P250-z5-multiseparation-smoke"


class Z5ChargedMultiseparationTests(unittest.TestCase):
    def test_exact_gate_covers_all_separations(self) -> None:
        gate = exact_mapping_gate()
        self.assertTrue(gate["passed"])
        self.assertEqual(gate["separations"], [1, 2, 3])
        self.assertEqual(gate["root_labels_checked"], 5850)
        self.assertEqual(gate["collapse_failures"], 0)
        self.assertEqual(gate["projection_failures"], 0)

    def test_hermitian_pair_is_exactly_real(self) -> None:
        origin = {1: complex(1, 2), 2: complex(3, -1), 3: complex(3, 1), 4: complex(1, -2)}
        xrow = {1: complex(-2, 1), 2: complex(2, 4), 3: complex(2, -4), 4: complex(-2, -1)}
        yrow = {1: complex(5, -3), 2: complex(-1, 2), 3: complex(-1, -2), 4: complex(5, 3)}
        self.assertEqual(hermitian_pair(origin, xrow, yrow, 1).imag, 0.0)
        self.assertEqual(hermitian_pair(origin, xrow, yrow, 2).imag, 0.0)

    def test_tiny_stream_retains_joint_72_coordinate_covariance(self) -> None:
        rows, analysis = run(4, 2, 1, 0.59274605079, 25011312220260901, 1, 0)
        self.assertEqual(len(rows), 2)
        self.assertEqual(len(FIELD_ORDER), 72)
        self.assertEqual(len(analysis["covariance_of_mean"]), 72)
        self.assertEqual(SMOKE_CAP, 5000)
        self.assertLess(analysis["pair_imaginary_max"], 1e-12)
        self.assertLess(analysis["DFT_conjugacy_max_abs"], 1e-12)

    def test_checked_smoke_applies_support_before_phase(self) -> None:
        payload = json.loads((RESULT / "response_4k.json").read_text())
        observed = score(payload, read_batches(RESULT / "response_4k.batches.csv"))
        checked = json.loads((RESULT / "score_4k.json").read_text())
        self.assertEqual(observed, checked)
        self.assertTrue(checked["separations"]["1"]["two_point_ready_abs_z_ge_2"])
        self.assertFalse(checked["separations"]["2"]["two_point_ready_abs_z_ge_2"])
        self.assertFalse(checked["separations"]["3"]["two_point_ready_abs_z_ge_2"])
        for row in checked["separations"].values():
            self.assertEqual(
                row["cubic_support_zero_score"]["decision_at_0.05"],
                "not_detected",
            )
            self.assertEqual(
                row["phase_closure"]["decision_status"],
                "not_interpretable_until_nonzero_support",
            )

    def test_existing_archive_audit_is_reproducible(self) -> None:
        observed = audit(
            ROOT / "results/huawei-20260830/P250-z5-charged-threepoint/score_1m.json",
            ROOT / "results/huawei-20260830/P250-z5-charged-threepoint/charged_threepoint_1m_merged.batches.csv",
            ROOT / "results/huawei-20260829/P226-norm5-chiral-fixedp/chiral_response.batches.csv",
        )
        checked = json.loads((RESULT / "existing_archive_reveal.json").read_text())
        self.assertEqual(observed, checked)
        self.assertFalse(checked["can_construct_same_operator_separation_ratio"])
        self.assertEqual(
            checked["existing_1m_cubic_support"]["decision_at_0.05"],
            "not_detected",
        )


if __name__ == "__main__":
    unittest.main()
