
from __future__ import annotations
import json
from pathlib import Path
import sys
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from pinson_arguin_kdv import (  # noqa: E402
    artifact_payload,
    primitive_k4_holomorphic_numeric,
    primitive_k4_holomorphic_series,
    primitive_numerator_extended,
)
from pinson_arguin_primitive import primitive_probability_direct  # noqa: E402
from score_p231_vacuum_kdv_sector import build_score  # noqa: E402


class PinsonArguinKdVTests(unittest.TestCase):
    def setUp(self) -> None:
        mp.mp.dps = 90

    def test_q_series_matches_numerical_wirtinger_derivative(self) -> None:
        tau = mp.mpc(mp.mpf("0.37"), mp.mpf("1.11"))
        for sector in ((1, 0), (0, 1), (1, 1), (2, -1)):
            series = primitive_k4_holomorphic_series(*sector, tau, dps=70)
            numeric = primitive_k4_holomorphic_numeric(*sector, tau, dps=70)
            self.assertLess(abs(series - numeric), mp.mpf("1e-58"))

    def test_q1_restricted_numerator_equals_probability(self) -> None:
        tau = mp.mpc(mp.mpf("0.41"), mp.mpf("0.93"))
        for sector in ((1, 0), (0, 1), (1, 1)):
            numerator = primitive_numerator_extended(
                *sector, tau, mp.conj(tau), dps=70
            )
            probability = primitive_probability_direct(*sector, tau, dps=70)
            self.assertLess(abs(numerator - probability), mp.mpf("1e-58"))

    def test_positive_rho_reflection_null_is_exact_numerically(self) -> None:
        payload = artifact_payload(dps=70)
        records = {record["id"]: record for record in payload["records"]}
        for identifier in ("pell_Dminus2_N30", "pell_Dplus1_N56"):
            record = records[identifier]
            response = [mp.mpf(value) for value in record["reflection_even_K4_plus_K4bar"]]
            self.assertLess(abs(response[1] - response[2]), mp.mpf("1e-58"))
            contrasts = [mp.mpf(value) for value in record["reflection_even_contrasts"]]
            self.assertLess(abs(contrasts[1]), mp.mpf("1e-58"))

    def test_committed_oracle_is_reproducible(self) -> None:
        path = ROOT / "predictions/p231_pinson_arguin_kdv_oracles_20260829.json"
        if not path.exists():
            self.skipTest("oracle is generated after the first focused test run")
        self.assertEqual(
            json.loads(path.read_text(encoding="utf-8")), artifact_payload(dps=80)
        )

    def test_retrospective_score_records_full_failure_and_C_compatibility(self) -> None:
        pilot = ROOT / "results/local-20260829/P156-square-bond-primitive-pilot/result.json"
        oracle = ROOT / "predictions/p231_pinson_arguin_kdv_oracles_20260829.json"
        if not oracle.exists():
            self.skipTest("oracle is generated after the first focused test run")
        score = build_score(pilot, oracle)
        self.assertEqual(score["decision"]["full_sector_vector"], "FAIL")
        self.assertEqual(score["decision"]["non_scalar_C_direction"], "COMPATIBLE")
        self.assertTrue(score["structural_predictions"]["C_same_sign_across_N30_N56"])
        self.assertAlmostEqual(
            score["structural_predictions"]["C_theory_ratio_N30_over_N56"],
            1.99068779989627,
            places=13,
        )
        self.assertTrue(
            all(
                abs(value) < 1e-50
                for value in score["structural_predictions"]["Q_reflection_null"]
            )
        )


if __name__ == "__main__":
    unittest.main()
