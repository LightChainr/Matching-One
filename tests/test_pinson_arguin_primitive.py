from __future__ import annotations

import json
from math import gcd
from pathlib import Path
import sys
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from pinson_arguin_primitive import (  # noqa: E402
    artifact_payload,
    canonical_primitive,
    engine_to_paper,
    paper_to_engine,
    primitive_probability_direct,
    primitive_probability_theta,
)


def assert_mpf_close(
    case: unittest.TestCase,
    value: mp.mpf,
    expected: str,
    tolerance: str = "1e-48",
) -> None:
    case.assertLess(abs(value - mp.mpf(expected)), mp.mpf(tolerance))


class PrimitiveConventionTests(unittest.TestCase):
    def test_engine_paper_sign_map_and_saturation(self) -> None:
        self.assertEqual(engine_to_paper((1, 0)), (1, 0))
        self.assertEqual(engine_to_paper((0, 1)), (0, 1))
        self.assertEqual(engine_to_paper((1, -1)), (1, 1))
        self.assertEqual(engine_to_paper((-2, 2)), (1, 1))
        self.assertEqual(paper_to_engine((1, 1)), (1, -1))

    def test_formula_rejects_zero_and_nonprimitive_labels(self) -> None:
        with self.assertRaisesRegex(ValueError, "zero winding"):
            canonical_primitive((0, 0))
        with self.assertRaisesRegex(ValueError, "coprime"):
            primitive_probability_direct(2, 0, 1j)


class PrimitiveProbabilityTests(unittest.TestCase):
    def setUp(self) -> None:
        mp.mp.dps = 100

    def test_direct_and_theta_formulas_are_independent_cross_checks(self) -> None:
        tau = mp.mpc(mp.mpf("0.37"), mp.mpf("1.11"))
        for sector in ((1, 0), (0, 1), (1, 1), (2, -1)):
            direct = primitive_probability_direct(*sector, tau, dps=80)
            theta = primitive_probability_theta(*sector, tau, dps=80)
            self.assertLess(abs(direct - theta), mp.mpf("1e-70"))

    def test_modular_T_and_S_laws(self) -> None:
        tau = mp.mpc(mp.mpf("0.37"), mp.mpf("1.11"))
        a, b = 2, -1
        original = primitive_probability_direct(a, b, tau, dps=80)
        transformed_T = primitive_probability_direct(
            a + b, b, tau + 1, dps=80
        )
        transformed_S = primitive_probability_direct(
            -b, a, -1 / tau, dps=80
        )
        self.assertLess(abs(original - transformed_T), mp.mpf("1e-70"))
        self.assertLess(abs(original - transformed_S), mp.mpf("1e-70"))

    def test_engine_basis_changes_match_paper_modular_laws(self) -> None:
        # P' = P U.  Repository windings transform as w' = U^-1 w.
        engine = (2, 1)
        paper = engine_to_paper(engine)

        # U_T=((1,1),(0,1)); U_T^-1 sends (u,v)->(u-v,v).
        engine_T = (engine[0] - engine[1], engine[1])
        self.assertEqual(engine_to_paper(engine_T), (paper[0] + paper[1], paper[1]))

        # U_S=((0,-1),(1,0)); U_S^-1 sends (u,v)->(v,-u).
        engine_S = (engine[1], -engine[0])
        self.assertEqual(engine_to_paper(engine_S), (-paper[1], paper[0]))

    def test_tau_i_oracles(self) -> None:
        tau = mp.mpc(0, 1)
        assert_mpf_close(
            self,
            primitive_probability_direct(1, 0, tau, dps=90),
            "0.16941543532134688938260796919875445000145337645375",
        )
        assert_mpf_close(
            self,
            primitive_probability_direct(1, 1, tau, dps=90),
            "0.020979928575590629661470611008187000992928553334531",
        )

    def test_half_sheared_sign_oracle(self) -> None:
        tau = mp.mpc(mp.mpf("0.5"), 1)
        assert_mpf_close(
            self,
            primitive_probability_direct(1, 0, tau, dps=90),
            "0.16815464971788045003554385835401755689282702879997",
        )
        equal = primitive_probability_direct(0, 1, tau, dps=90)
        assert_mpf_close(
            self,
            equal,
            "0.10005678718797952632971693447464511165342474778502",
        )
        self.assertLess(
            abs(equal - primitive_probability_direct(1, 1, tau, dps=90)),
            mp.mpf("1e-80"),
        )
        assert_mpf_close(
            self,
            primitive_probability_direct(1, -1, tau, dps=90),
            "0.0015190922810096831710181718634974911858075569375806",
        )

    def test_hexagonal_three_sector_oracle(self) -> None:
        tau = mp.mpc(mp.mpf("0.5"), mp.sqrt(3) / 2)
        expected = "0.12166379946598032273800506616100306406127073342455"
        values = [
            primitive_probability_direct(*sector, tau, dps=90)
            for sector in ((1, 0), (0, 1), (1, 1))
        ]
        for value in values:
            assert_mpf_close(self, value, expected)
        assert_mpf_close(
            self,
            primitive_probability_direct(1, -1, tau, dps=90),
            "0.00096713702243181880247990504669184802012174757920661",
        )

    def test_small_pell_oracles(self) -> None:
        rows = (
            (
                mp.mpc(mp.mpf("0.5"), mp.mpf(5) / 6),
                "0.11072917769038501312751860913865363519163224173796",
                "0.1272155037499346429480553869141822763171677841226",
            ),
            (
                mp.mpc(mp.mpf("0.5"), mp.mpf(7) / 8),
                "0.12470059657156361356306944379819518283006971050975",
                "0.12015211579226502308743616256051247672799257647673",
            ),
        )
        for tau, horizontal, diagonal in rows:
            assert_mpf_close(
                self,
                primitive_probability_direct(1, 0, tau, dps=90),
                horizontal,
            )
            for sector in ((0, 1), (1, 1)):
                assert_mpf_close(
                    self,
                    primitive_probability_direct(*sector, tau, dps=90),
                    diagonal,
                )

    def test_square_normalization_and_cross_duality(self) -> None:
        tau = mp.mpc(0, 1)
        primitive_sum = mp.mpf(0)
        cutoff = 10
        for a in range(0, cutoff + 1):
            for b in range(-cutoff, cutoff + 1):
                if a == 0 and b <= 0:
                    continue
                if gcd(a, abs(b)) != 1:
                    continue
                primitive_sum += primitive_probability_direct(
                    a, b, tau, dps=70
                )
        assert_mpf_close(
            self,
            primitive_sum,
            "0.38094744914033735446061273329420244691096663128102",
            "1e-45",
        )
        cross = (1 - primitive_sum) / 2
        assert_mpf_close(
            self,
            cross,
            "0.30952627542983132276969363335289877654451668435949",
            "1e-45",
        )


class PrimitiveArtifactTests(unittest.TestCase):
    def test_machine_readable_artifact_matches_evaluator(self) -> None:
        artifact = json.loads(
            (
                ROOT
                / "predictions"
                / "p156_pinson_arguin_baselines_20260829.json"
            ).read_text(encoding="utf-8")
        )
        generated = artifact_payload(dps=90)
        self.assertEqual(artifact, generated)
        self.assertEqual(artifact["issue"], 156)
        self.assertIn("rank-1", artifact["convention"]["event"])


if __name__ == "__main__":
    unittest.main()
