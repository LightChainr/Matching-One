import json
from pathlib import Path
import unittest

from p334_tm_configuration_cross_switch import build_result


ROOT = Path(__file__).resolve().parents[1]


class P334TMConfigurationCrossSwitchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result()

    def test_strict_mark_switch_is_refuted_minimally(self):
        self.assertEqual(
            self.result["universal_strict_injection_status"], "refuted"
        )
        gate = self.result["minimal_counterexample_gate"]
        self.assertEqual(gate["hard_rows_below_N6"], 0)
        self.assertEqual(gate["minimal_N"], 6)
        self.assertEqual(gate["minimal_rows"], 4)
        for row in self.result["rows"]:
            self.assertEqual(row["forced_alignment_positive_targets"], 0)
            self.assertEqual(
                row["operation_audits"]["mark_only"]["Hall_deficiency"],
                1032,
            )

    def test_phase_and_replica_are_exact_lost_labels(self):
        for row in self.result["rows"]:
            phase = row["phase_information_loss"]
            self.assertEqual(phase["bare_target"]["maximum_fiber"], 24)
            self.assertEqual(
                phase["target_plus_translation_phase"]["maximum_fiber"], 4
            )
            self.assertEqual(
                phase["target_plus_phase_and_replica"]["maximum_fiber"], 1
            )

    def test_union_and_one_carrier_moves_still_fail(self):
        union = set()
        one_carrier = set()
        for row in self.result["rows"]:
            audits = row["operation_audits"]
            union.add(audits["union_preserving"]["Hall_deficiency"])
            one_carrier.add(
                audits["one_carrier_transport"]["Hall_deficiency"]
            )
        self.assertEqual(union, {840, 852})
        self.assertEqual(one_carrier, {276, 336})

    def test_two_carrier_transport_saturates_minimal_rows(self):
        for row in self.result["rows"]:
            audit = row["operation_audits"]["two_carrier_transport"]
            self.assertEqual(audit["source_tokens"], 1152)
            self.assertEqual(audit["maximum_matching"], 1152)
            self.assertEqual(audit["Hall_deficiency"], 0)
            self.assertTrue(audit["collision_free_injection_exists"])

    def test_checked_artifact_reproduces(self):
        checked = json.loads(
            (
                ROOT
                / "results/p334-tm-configuration-cross-switch/latest.json"
            ).read_text()
        )
        self.assertEqual(checked, self.result)


if __name__ == "__main__":
    unittest.main()
