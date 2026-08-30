from __future__ import annotations

import copy
import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import pc_confidence_preregistration as prereg  # noqa: E402


class PcConfidencePreregistrationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.plan = prereg.build_plan()
        self.records = prereg.synthetic_fixture(self.plan)

    def test_frozen_exact_constants(self) -> None:
        contract = self.plan["statistical_contract"]
        self.assertEqual(contract["familywise_alpha"], "1/1000000")
        self.assertEqual(contract["per_run_alpha"], "1/6000000")
        self.assertEqual(contract["minimum_successes"], 373)
        self.assertEqual(contract["trials_per_attempt"], 400)

    def test_synthetic_fixture_passes_and_is_not_empirical(self) -> None:
        audit = prereg.audit_records(self.plan, self.records)
        self.assertEqual(audit["record_count"], 6)
        self.assertEqual(audit["unique_side_attempts"], 6)
        self.assertEqual(audit["accepted_records"], 2)
        self.assertTrue(all(row["tested_parameter"] == "synthetic-only" for row in self.records))

    def test_plan_tampering_invalidates_existing_records(self) -> None:
        changed = copy.deepcopy(self.plan)
        changed["statistical_contract"]["trials_per_attempt"] = 401
        with self.assertRaisesRegex(ValueError, "plan digest mismatch"):
            prereg.audit_records(changed, self.records)

    def test_duplicate_side_attempt_fails(self) -> None:
        duplicate = copy.deepcopy(self.records[0])
        duplicate["record_id"] = "different-id"
        duplicate["stream_domain"] = "different-domain"
        duplicate["data_digest"] = prereg.sha256_text("different-data")
        with self.assertRaisesRegex(ValueError, "duplicate side/attempt"):
            prereg.audit_records(self.plan, self.records + [duplicate])

    def test_duplicate_domain_and_data_fail(self) -> None:
        for field, message in (("stream_domain", "duplicate stream_domain"), ("data_digest", "duplicate data_digest")):
            changed = copy.deepcopy(self.records)
            changed[1][field] = changed[0][field]
            with self.assertRaisesRegex(ValueError, message):
                prereg.audit_records(self.plan, changed)

    def test_exploration_digest_reuse_fails_after_redigesting_plan(self) -> None:
        changed_plan = copy.deepcopy(self.plan)
        changed_plan["forbidden_exploration_data_digests"] = [self.records[0]["data_digest"]]
        changed_records = copy.deepcopy(self.records)
        digest = prereg.plan_digest(changed_plan)
        for record in changed_records:
            record["plan_digest"] = digest
        with self.assertRaisesRegex(ValueError, "reuses an exploration digest"):
            prereg.audit_records(changed_plan, changed_records)

    def test_graph_counts_and_phase_fail_closed(self) -> None:
        mutations = [
            ("graph", "matching", "graph/side mismatch"),
            ("trials", 399, "trial count differs"),
            ("successes", 401, "invalid success count"),
            ("phase", "exploration", "not a final trial"),
        ]
        for field, value, message in mutations:
            changed = copy.deepcopy(self.records)
            changed[0][field] = value
            with self.assertRaisesRegex(ValueError, message):
                prereg.audit_records(self.plan, changed)

    def test_cutoff_matches_exact_tail_on_both_sides(self) -> None:
        audit = prereg.audit_records(self.plan, self.records)
        by_successes = {row["successes"]: row["accepted"] for row in audit["decisions"]}
        self.assertFalse(by_successes[372])
        self.assertTrue(by_successes[373])

    def test_checked_in_results_reproduce(self) -> None:
        artifact = prereg.build_artifact()
        checked_json = json.loads(
            (ROOT / "results/pc-confidence-preregistration/latest.json").read_text(encoding="utf-8")
        )
        checked_markdown = (
            ROOT / "results/pc-confidence-preregistration/latest.md"
        ).read_text(encoding="utf-8")
        self.assertEqual(checked_json, artifact)
        self.assertEqual(checked_markdown, prereg.render_markdown(artifact))


if __name__ == "__main__":
    unittest.main()
