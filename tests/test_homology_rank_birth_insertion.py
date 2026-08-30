from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import homology_rank_birth_insertion as oracle  # noqa: E402


class HomologyRankBirthInsertionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.artifact = oracle.build_artifact()
        cls.rows = {row["id"]: row for row in cls.artifact["geometries"]}

    def test_pointwise_gates_include_direct_rank_two_jump(self) -> None:
        geometry = oracle.axis_integer_torus(2)
        # The exhaustive oracle freezes this as its first 0->2 representative.
        active = oracle._mask_active(6, geometry.n)
        insertion = oracle.rank_birth_insertion(geometry, active, 0)
        self.assertEqual((insertion["rank_before"], insertion["rank_after"]), (0, 2))
        self.assertEqual(insertion["delta_rank"], 2)
        self.assertEqual(insertion["gate_0_to_1"], 1)
        self.assertEqual(insertion["gate_1_to_2"], 1)
        self.assertEqual([birth["ell"] for birth in insertion["births"]], [None, None])

    def test_exact_russo_polynomials_close_on_every_geometry(self) -> None:
        for row in self.rows.values():
            coefficients = row["power_basis_coefficients_ascending"]
            self.assertEqual(
                coefficients["M_prime_direct"], coefficients["f_01_plus_f_12"]
            )
        self.assertEqual(
            self.rows["axis-L4-fixed-root"]["p_equals_half"],
            {"M_prime": "4209/1024", "f_01": "5297/2048", "f_12": "3121/2048"},
        )

    def test_endpoint_line_and_index_semantics(self) -> None:
        representatives = self.rows["gaussian-2-1"]["representatives"]
        first = representatives["0->1"]["insertion"]["births"][0]
        second = representatives["1->2"]["insertion"]["births"][0]
        self.assertEqual(first["ell_role"], "new_rank_one_image_after_0_to_1")
        self.assertEqual(second["ell_role"], "rank_one_plateau_line_before_1_to_2")
        self.assertIsNotNone(first["ell"])
        self.assertIsNotNone(second["ell"])
        self.assertEqual(first["iota"], 1)
        self.assertEqual(second["iota"], 1)
        for row in self.rows.values():
            for counts in row["line_index_counts"].values():
                self.assertFalse(any("iota=" in key and "iota=1" not in key for key in counts))

    def test_homology_h4_is_exact_and_local_mark_is_attached(self) -> None:
        geometry = oracle.gaussian_integer_torus(2, 1)
        mark = oracle.homology_h4_mark(geometry, (0, 1))
        self.assertEqual(mark, {"physical_vector": "-1,2", "cos4": "-7/25", "sin4": "24/25"})
        local = self.rows["axis-L4-fixed-root"]["local_landing_h4_raw_sums"]
        self.assertEqual(local["0_to_1"]["h4"], 8256)
        self.assertEqual(local["1_to_2"]["h4"], 6688)
        self.assertGreater(local["0_to_1"]["landed"], 0)
        self.assertGreater(local["1_to_2"]["landed"], 0)

    def test_axis_l4_has_simultaneous_births_and_landing_is_not_rank_birth(self) -> None:
        row = self.rows["axis-L4-fixed-root"]
        self.assertEqual(row["transition_counts"]["0->2"], 4624)
        first_gate_total = sum(row["birth_bernstein_counts_by_other_occupancy"]["0_to_1"])
        second_gate_total = sum(row["birth_bernstein_counts_by_other_occupancy"]["1_to_2"])
        self.assertEqual(first_gate_total, 84752)
        self.assertEqual(second_gate_total, 49936)
        self.assertLess(row["local_landing_h4_raw_sums"]["0_to_1"]["landed"], first_gate_total)
        self.assertLess(row["local_landing_h4_raw_sums"]["1_to_2"]["landed"], second_gate_total)

    def test_checked_in_artifacts_reproduce(self) -> None:
        directory = ROOT / "results/homology-rank-birth-insertion"
        checked_json = json.loads((directory / "latest.json").read_text(encoding="utf-8"))
        checked_markdown = (directory / "latest.md").read_text(encoding="utf-8")
        self.assertEqual(checked_json, self.artifact)
        self.assertEqual(checked_markdown, oracle.render_markdown(self.artifact))


if __name__ == "__main__":
    unittest.main()

