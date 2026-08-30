from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import p334_orbit_flux_phase_diagram as phase  # noqa: E402
from p334_n13_multiorbit_flux import P_REF, _evaluate_incidence, exact_census  # noqa: E402


class P334OrbitFluxPhaseDiagramTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = phase.build_certificate()

    def test_monomial_expansion_matches_bernstein_evaluation(self) -> None:
        for a, b in ((3, 2), (4, 1)):
            census = exact_census(a, b, include_direct_rank2=True)
            for label in phase.LABELS:
                for source, key in (("birth", "birth_edges"), ("exit", "exit_edges")):
                    poly = phase._bernstein_edge_poly(census, label, key)
                    self.assertEqual(
                        phase._eval(poly, P_REF),
                        _evaluate_incidence(census, label, source, P_REF),
                    )

    def test_all_interior_roots_are_isolated_and_complete(self) -> None:
        for geometry in self.payload["geometries"].values():
            for name in (
                "axis_birth",
                "axis_exit",
                "diagonal_birth",
                "diagonal_exit",
                "birth_character_total",
                "exit_character_total",
            ):
                self.assertEqual(geometry["roots"][name]["interior_root_count"], 0)
            for name in ("axis_net", "diagonal_net", "total_net"):
                root_set = geometry["roots"][name]
                self.assertEqual(root_set["interior_root_count"], 1)
                row = root_set["roots"][0]
                self.assertLessEqual(
                    Fraction(row.get("width", "0")), Fraction(1, 1 << 111)
                )

    def test_root_order_reverses_but_p_ref_stays_in_reinforcing_window(self) -> None:
        n13 = self.payload["geometries"]["N13"]["signed_share_singularities"]
        n17 = self.payload["geometries"]["N17"]["signed_share_singularities"]
        root = lambda rows: phase._root_midpoint(rows[0])
        self.assertLess(root(n13["diagonal_share_zero"]), P_REF)
        self.assertLess(P_REF, root(n13["axis_share_zero"]))
        self.assertLess(root(n17["axis_share_zero"]), P_REF)
        self.assertLess(P_REF, root(n17["diagonal_share_zero"]))
        for geometry in self.payload["geometries"].values():
            interval = next(
                row
                for row in geometry["phase_intervals"]
                if Fraction(row["lower"]) < P_REF < Fraction(row["upper"])
            )
            self.assertEqual(interval["orbit_contributions"], "reinforce")

    def test_point_stability_is_not_slope_stability(self) -> None:
        n13 = self.payload["geometries"]["N13"]["p_ref_metrics"]
        n17 = self.payload["geometries"]["N17"]["p_ref_metrics"]
        self.assertAlmostEqual(
            float(n13["axis_signed_share"]), float(n17["axis_signed_share"]), delta=0.01
        )
        self.assertLess(float(n13["axis_signed_share_slope_dp"]), 0)
        self.assertGreater(float(n17["axis_signed_share_slope_dp"]), 0)
        for point in (n13, n17):
            self.assertLess(float(point["activity_cancellation"]["axis"]), 0.08)
            self.assertLess(float(point["activity_cancellation"]["diagonal"]), 0.08)

    def test_checked_in_artifacts_reproduce(self) -> None:
        expected_json = json.loads(
            (ROOT / "results/p334-orbit-flux-phase-diagram/latest.json").read_text(
                encoding="utf-8"
            )
        )
        expected_markdown = (
            ROOT / "results/p334-orbit-flux-phase-diagram/latest.md"
        ).read_text(encoding="utf-8")
        self.assertEqual(expected_json, self.payload)
        self.assertEqual(expected_markdown, phase.render_markdown(self.payload) + "\n")


if __name__ == "__main__":
    unittest.main()
