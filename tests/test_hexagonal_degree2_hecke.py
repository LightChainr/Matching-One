import json
import sys
import unittest
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from derive_hexagonal_degree2_hecke import (  # noqa: E402
    eisenstein_qseries,
    fingerprint,
)


class HexagonalDegree2HeckeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = fingerprint(90)
        cls.artifact = json.loads(
            (ROOT / "predictions" / "hexagonal_degree2_hecke_20260829.json")
            .read_text(encoding="utf-8")
        )

    def test_direct_qseries_fingerprint(self) -> None:
        errors = self.payload["independent_checks"]["fingerprint_errors"]
        self.assertLess(mp.mpf(errors["E4_parent_zero_abs"]), mp.mpf("1e-80"))
        self.assertLess(mp.mpf(errors["E4_child_phase_max_abs"]), mp.mpf("1e-68"))
        self.assertLess(mp.mpf(errors["E6_child_ratio_max_abs"]), mp.mpf("1e-68"))

    def test_modular_and_hecke_cross_checks(self) -> None:
        for weight in (4, 6):
            checks = self.payload["independent_checks"][f"weight_{weight}"]
            self.assertLess(mp.mpf(checks["S_then_T_error_abs"]), mp.mpf("1e-80"))
            self.assertLess(mp.mpf(checks["reflection_error_abs"]), mp.mpf("1e-80"))
            self.assertLess(mp.mpf(checks["Hecke_residual_abs"]), mp.mpf("1e-80"))

    def test_raw_e6_ratios(self) -> None:
        rows = self.payload["numerical_qseries"]["weight_6"]
        parent = mp.mpc(rows["parent"]["raw_Ek"]["real"], rows["parent"]["raw_Ek"]["imag"])
        targets = [mp.mpf(11) / 32, mp.mpf(22), mp.mpf(22)]
        for row, target in zip(rows["children"], targets):
            child = mp.mpc(row["raw_Ek"]["real"], row["raw_Ek"]["imag"])
            self.assertLess(abs(child / parent - target), mp.mpf("1e-68"))

    def test_qseries_rejects_other_weights(self) -> None:
        with self.assertRaises(ValueError):
            eisenstein_qseries(8, mp.j)

    def test_frozen_artifact_is_reproducible(self) -> None:
        self.assertEqual(self.artifact, self.payload)


if __name__ == "__main__":
    unittest.main()
