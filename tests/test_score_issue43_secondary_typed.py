
from __future__ import annotations
from pathlib import Path
import shutil
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_issue43_secondary_typed as typed  # noqa: E402


def payload(last_status: str = "READY_AWAITING_DERIVATIVE_TARGET") -> dict:
    statuses = typed.FIXED_STATUSES + [last_status]
    return {
        "protocol": "Issue #43 frozen secondary scoring ledger",
        "status": "fixed secondary scores only; no target refit",
        "sizes": [185, 265],
        "stage_order": list(typed.STAGE_ORDER),
        "stages": [
            {"order": index + 1, "name": name, "status": status}
            for index, (name, status) in enumerate(zip(typed.STAGE_ORDER, statuses))
        ],
        "excluded_models": [{
            "name": "V_<1,3>_N^-4/3",
            "status": "EXCLUDED_INVALIDATED_WRONG_KAC_BRANCH",
            "scored": False,
        }],
        "raw_data_boundary": "frozen",
    }


class Issue43SecondaryTypedTests(unittest.TestCase):
    def copied_root(self) -> tuple[tempfile.TemporaryDirectory, Path]:
        directory = tempfile.TemporaryDirectory()
        root = Path(directory.name)
        destination = root / typed.SEMANTIC_GATE
        destination.parent.mkdir(parents=True)
        shutil.copy(ROOT / typed.SEMANTIC_GATE, destination)
        return directory, root

    def test_gate_registers_three_exact_identities(self) -> None:
        _, validated = typed.load_semantic_gate(ROOT)
        self.assertEqual(list(validated), ["DeltaM", "DeltaS", "P4_S_prime"])
        self.assertTrue(all(
            (item["transform"].scale, item["transform"].offset) == (1.0, 0.0)
            for item in validated.values()
        ))

    def test_payload_is_unchanged_before_annotation(self) -> None:
        frozen = payload()
        result = typed.score_typed(
            ROOT, {}, Path("x17"), Path("p48"), runner=lambda *_: frozen
        )
        semantics = result.pop("observable_semantics")
        self.assertIs(result, frozen)
        self.assertEqual(result, payload())
        self.assertEqual(semantics["stage_order"], typed.STAGE_ORDER)

    def test_not_scorable_stage_cannot_be_promoted(self) -> None:
        bad = payload()
        bad["stages"][3]["status"] = "SCORED_FROZEN_NO_REFIT"
        with self.assertRaisesRegex(ValueError, "stage status"):
            typed.score_typed(ROOT, {}, Path("x17"), Path("p48"), runner=lambda *_: bad)

    def test_excluded_model_cannot_be_scored(self) -> None:
        bad = payload()
        bad["excluded_models"][0]["scored"] = True
        with self.assertRaisesRegex(ValueError, "excluded-model"):
            typed.score_typed(ROOT, {}, Path("x17"), Path("p48"), runner=lambda *_: bad)


if __name__ == "__main__":
    unittest.main()
