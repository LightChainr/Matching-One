import copy
import importlib.util
import json
import pathlib
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "score_matching_odd_synthesis", ROOT / "scripts" / "score_matching_odd_synthesis.py"
)
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


class MatchingOddSynthesisTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ledger = json.loads((ROOT / "results" / "evidence-ledger" / "latest.json").read_text())

    def test_committed_ledger_reproduces_issue212_numbers(self):
        result = MOD.synthesize(self.ledger)
        zero = result["joint_scores"]["zero_effect"]
        h4 = result["joint_scores"]["fixed_H4"]
        self.assertAlmostEqual(zero["chi_square"], 31.18573555150965)
        self.assertEqual(zero["dof"], 4)
        self.assertAlmostEqual(zero["chi_square_survival_p"], 2.805595267905808e-6)
        self.assertAlmostEqual(h4["chi_square"], 3.4622795373044296)
        self.assertEqual(h4["dof"], 4)
        self.assertAlmostEqual(h4["chi_square_survival_p"], 0.48363695393249573)
        self.assertAlmostEqual(h4["nlpd"], -35.7946059274312)
        self.assertAlmostEqual(zero["nlpd"], -21.988187137702105)
        self.assertAlmostEqual(
            result["predictive_comparison"]["delta_nlpd_fixed_H4_minus_zero_effect"],
            -13.806418789729096,
        )
        self.assertFalse(result["governance"]["adds_new_primary_evidence"])

    def test_nonprimary_or_wrong_channel_is_rejected(self):
        for mutation, message in (({"role": "sensitivity"}, "not primary"), ({"channel": {"source": "matching_odd", "target": "other"}}, "matching_odd")):
            with self.subTest(mutation=mutation):
                ledger = copy.deepcopy(self.ledger)
                block = next(row for row in ledger["blocks"] if row["id"] == MOD.BLOCK_SPECS[0]["id"])
                block.update(mutation)
                with self.assertRaisesRegex(ValueError, message):
                    MOD.synthesize(ledger)

    def test_raw_data_group_contract_and_score_registration_are_enforced(self):
        ledger = copy.deepcopy(self.ledger)
        block = next(row for row in ledger["blocks"] if row["id"] == MOD.BLOCK_SPECS[1]["id"])
        block["raw_data_group"] = MOD.BLOCK_SPECS[0]["raw_data_group"]
        with self.assertRaisesRegex(ValueError, "raw_data_group"):
            MOD.synthesize(ledger)

        ledger = copy.deepcopy(self.ledger)
        block = next(row for row in ledger["blocks"] if row["id"] == MOD.BLOCK_SPECS[1]["id"])
        del block["scores"]["H4_norm5"]
        with self.assertRaisesRegex(ValueError, "lacks registered score"):
            MOD.synthesize(ledger)


if __name__ == "__main__":
    unittest.main()
