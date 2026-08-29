from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p275_double_projector_selection import (  # noqa: E402
    A_TOP,
    ALEXANDER,
    EVEN_TANGENT,
    I3,
    ODD_TANGENT,
    artifact_payload,
    dot,
    jordan_parity_certificate,
    matmul,
    matvec,
    scale,
    thermal_q4_certificate,
    vacuum_kdv_certificate,
)


class P275DoubleProjectorSelectionTests(unittest.TestCase):
    def test_alexander_involution_and_rank_lines(self) -> None:
        self.assertEqual(matmul(ALEXANDER, ALEXANDER), I3)
        self.assertEqual(matvec(ALEXANDER, EVEN_TANGENT), EVEN_TANGENT)
        self.assertEqual(matvec(ALEXANDER, ODD_TANGENT), scale(-1, ODD_TANGENT))
        self.assertEqual(sum(EVEN_TANGENT), 0)
        self.assertEqual(sum(ODD_TANGENT), 0)
        self.assertEqual(dot(A_TOP, EVEN_TANGENT), 0)
        self.assertEqual(dot(A_TOP, ODD_TANGENT), 2)

    def test_vacuum_kdv_is_exactly_annihilated(self) -> None:
        certificate = vacuum_kdv_certificate()
        self.assertEqual(certificate["status"], "EXACT_ZERO")
        self.assertEqual(certificate["response_direction"], ["1", "-2", "1"])
        self.assertEqual(certificate["identities"]["A_top_response"], "0")

    def test_thermal_q4_is_the_odd_rank_line(self) -> None:
        certificate = thermal_q4_certificate()
        self.assertEqual(certificate["status"], "CONDITIONAL_PARITY_THEOREM")
        self.assertEqual(certificate["response_direction"], ["-1", "0", "1"])
        self.assertEqual(certificate["identities"]["delta_P1"], "0")
        self.assertEqual(certificate["identities"]["delta_P2_plus_delta_P0"], "0")

    def test_jordan_partner_cannot_mix_with_bottom_under_involution(self) -> None:
        certificate = jordan_parity_certificate()
        self.assertEqual(certificate["derived_involution"], [["-1", "0"], ["0", "-1"]])
        self.assertEqual(certificate["status"], "EXACT_GIVEN_GRADING")

    def test_committed_certificate_is_reproducible(self) -> None:
        path = ROOT / "results/exact-double-projector-selection/latest.json"
        self.assertEqual(json.loads(path.read_text(encoding="utf-8")), artifact_payload())


if __name__ == "__main__":
    unittest.main()

