from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import rank_birth_parity_channels as oracle  # noqa: E402


class RankBirthParityChannelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.artifact = oracle.build_artifact()
        cls.rows = {row["id"]: row for row in cls.artifact["geometries"]}

    def test_complement_reversal_exchanges_the_two_gates(self) -> None:
        geometry = oracle.geometry_specs()[1]["geometry"]
        active = oracle._mask_active(6, geometry.n)
        primal, matching, _ = oracle.complement_pair(geometry, active, 0)
        self.assertEqual((primal["rank_before"], primal["rank_after"]), (0, 1))
        self.assertEqual((matching["rank_before"], matching["rank_after"]), (1, 2))
        self.assertEqual(
            (primal["gate_01"], primal["gate_12"]),
            (matching["gate_12"], matching["gate_01"]),
        )
        self.assertEqual(primal["ell"], matching["ell"])
        self.assertEqual(primal["iota"], matching["iota"])

    def test_all_environment_and_coefficient_residuals_are_zero(self) -> None:
        for row in self.rows.values():
            self.assertEqual(row["complement_failure_counts"], {})
            self.assertFalse(any(row["coefficient_parity_failures"].values()))

    def test_even_is_matching_derivative_and_odd_is_lifetime_derivative(self) -> None:
        axis = self.rows["axis-L4-fixed-root"]
        self.assertEqual(axis["p_equals_half_primal"]["even"], "4209/1024")
        self.assertEqual(axis["p_equals_half_primal"]["odd"], "-17/16")
        self.assertEqual(axis["exact_identifications"]["even"], "f_01+f_12=M_prime")
        self.assertIn("-partial_p Prob(rank=1)", axis["exact_identifications"]["odd"])

    def test_simultaneous_birth_is_even_and_has_no_odd_ambiguity(self) -> None:
        axis = self.rows["axis-L4-fixed-root"]
        self.assertEqual(axis["p_equals_half_primal"]["simultaneous_0_to_2"], "289/2048")
        strict_line_events = sum(axis["line_index_pair_counts"].values())
        paired = axis["paired_transition_counts"]
        self.assertEqual(
            strict_line_events,
            paired["0->1 | 1->2"] + paired["1->2 | 0->1"],
        )

    def test_projective_spin4_character_has_exact_gaussian_sign(self) -> None:
        for x, y in ((1, 0), (1, 1), (2, 1), (3, -2)):
            character = oracle.spin4_character((x, y))
            self.assertEqual(oracle.spin4_character((-x, -y)), character)
            self.assertEqual(oracle.spin4_character((-y, x)), character)
            self.assertEqual(
                oracle.spin4_character((x - y, x + y)),
                (-character[0], -character[1]),
            )

    def test_line_and_local_h4_odd_channels_are_nonzero(self) -> None:
        axis = self.rows["axis-L4-fixed-root"]["p_equals_half_primal"]
        self.assertEqual(axis["line_cos4_odd"], "-109/128")
        self.assertEqual(axis["local_h4_odd"], "-49/1024")
        gaussian = self.rows["gaussian-2-1"]["p_equals_half_primal"]
        self.assertEqual(gaussian["line_sin4_odd"], "-3/5")

    def test_checked_in_artifacts_reproduce(self) -> None:
        directory = ROOT / "results/rank-birth-parity-channels"
        checked_json = json.loads((directory / "latest.json").read_text(encoding="utf-8"))
        checked_markdown = (directory / "latest.md").read_text(encoding="utf-8")
        self.assertEqual(checked_json, self.artifact)
        self.assertEqual(checked_markdown, oracle.render_markdown(self.artifact))


if __name__ == "__main__":
    unittest.main()

