from __future__ import annotations

import importlib.util
import json
from fractions import Fraction
from pathlib import Path
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "p220_q4_logtorus_eta_cocycle.py"
SPEC = importlib.util.spec_from_file_location("p220_eta", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class Q4LogTorusEtaCocycleTests(unittest.TestCase):
    def test_fixed_q4_path_recovers_existing_Ward_coefficient(self) -> None:
        h = Fraction(5, 8)
        self.assertEqual(MODULE.fixed_q4_coefficient(h), Fraction(493, 96))
        self.assertEqual(MODULE.fixed_q4_coefficient_derivative(h), Fraction(3637, 540))
        self.assertEqual(
            MODULE.fixed_q4_coefficient_derivative(h) / MODULE.fixed_q4_coefficient(h),
            Fraction(29096, 22185),
        )

    def test_exact_eta_CM_ratios(self) -> None:
        with mp.workdps(80):
            eta_i = MODULE.eta(mp.j, dps=80)
            eta_2i = MODULE.eta(2 * mp.j, dps=80)
            eta_shear = MODULE.eta((1 + mp.j) / 2, dps=80)
            self.assertLess(abs(eta_2i / eta_i - mp.power(2, -mp.mpf(3) / 8)), mp.mpf("1e-70"))
            self.assertLess(abs(abs(eta_shear / eta_i) - mp.power(2, mp.mpf(1) / 4)), mp.mpf("1e-70"))

    def test_frozen_cocycle_ratio_is_minus_three_halves(self) -> None:
        with mp.workdps(80):
            eta_i = MODULE.eta(mp.j, dps=80)
            left = 2 * mp.log(abs(MODULE.eta(2 * mp.j, dps=80) / eta_i))
            right = 2 * mp.log(abs(MODULE.eta((1 + mp.j) / 2, dps=80) / eta_i))
            self.assertLess(abs(left / right + mp.mpf(3) / 2), mp.mpf("1e-70"))

    def test_hexagonal_bottom_and_energy_derivative_share_E4_zero(self) -> None:
        with mp.workdps(80):
            rho = (1 + mp.sqrt(3) * mp.j) / 2
            self.assertLess(abs(MODULE.e4(rho, dps=80)), mp.mpf("1e-70"))

    def test_artifact_is_reproducible(self) -> None:
        expected = json.loads((ROOT / "predictions" / "p220_q4_logtorus_eta_cocycle_20260829.json").read_text())
        self.assertEqual(MODULE.analyze(90), expected)


if __name__ == "__main__":
    unittest.main()
