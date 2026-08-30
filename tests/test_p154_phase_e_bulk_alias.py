from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p154_phase_e_bulk_alias import METRICS, SIZE_ORDER  # noqa: E402


class P154PhaseEBulkAliasTests(unittest.TestCase):
    def test_contract_covers_both_complete_lineages(self):
        self.assertEqual(SIZE_ORDER, (65, 130, 260, 520, 85, 170, 340, 680))
        self.assertEqual(
            METRICS,
            ("J_top_even_birth_U", "J_bulk_binomial_energy_U"),
        )

    def test_committed_alias_certificate_is_machine_zero(self):
        result = json.loads(
            (ROOT / "results/p154-phase-e-bulk-energy-alias/latest.json").read_text()
        )
        certificate = result["alias_certificate"]
        self.assertEqual(
            certificate["decision"],
            "exact_coordinate_alias_not_an_independent_bulk_direction",
        )
        self.assertLess(certificate["maximum_full_point_absolute_residual"], 1e-10)
        self.assertLess(
            certificate["maximum_explicit_delete_one_absolute_residual"], 1e-8
        )
        self.assertLess(certificate["maximum_difference_from_committed_U"], 1e-10)

    def test_no_independent_local_mark_is_claimed(self):
        result = json.loads(
            (ROOT / "results/p154-phase-e-bulk-energy-alias/latest.json").read_text()
        )
        self.assertFalse(result["schema_audit"]["independent_local_singlet_mark_present"])
        self.assertEqual(result["schema_audit"]["unexpected_fields"], {})


if __name__ == "__main__":
    unittest.main()
