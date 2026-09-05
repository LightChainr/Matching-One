
from __future__ import annotations
import copy
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_p50_fullcurve_n290_typed as typed  # noqa: E402


def frozen_result(contract: dict) -> dict:
    observation = {name: float(index) for index, name in enumerate(contract["feature_order"])}
    observation["p0"] = 0.5
    return {
        "schema": typed.FROZEN_SCHEMA,
        "status": typed.FROZEN_STATUS,
        "scoring_order": contract["scoring_order"],
        "observations": {"N145": observation.copy(), "N290": observation.copy()},
        "primary_deltaM_transfer": {"levels": [0.0, 0.025, 0.05]},
        "slope": {"raw_asymptotic_baseline": {}},
        "root": {"raw_minus_one_quarter_baseline": {}},
        "p48_diagnostics": [{"metric": name} for name in typed.P4_ORDER],
        "covariance_rule": typed.FROZEN_COVARIANCE_RULE,
        "provenance": {"fixture": True},
    }


class P50FullcurveTypedTests(unittest.TestCase):
    def paths(self) -> tuple[Path, Path, Path, Path, Path]:
        return (
            Path("parent.csv"), Path("parent.json"), Path("child.csv"),
            Path("child.json"),
            ROOT / "predictions/p49_slope_two_sector_145_290_20260828.yaml",
        )

    def test_contract_registers_four_identity_maps(self) -> None:
        contract, validated = typed.load_contract(ROOT)
        self.assertEqual(list(validated), ["thermal_even", "matching_function", "P4_S", "P4_D"])
        self.assertEqual(contract["sizes_in_order"], [145, 290])
        self.assertTrue(all(
            (row["transform"].scale, row["transform"].offset) == (1.0, 0.0)
            for row in validated.values()
        ))

    def test_delegated_result_is_numerically_unchanged(self) -> None:
        contract, _ = typed.load_contract(ROOT)
        frozen = frozen_result(contract)
        result = typed.score_typed(ROOT, *self.paths(), runner=lambda *_: copy.deepcopy(frozen))
        semantics = result.pop("observable_semantics")
        self.assertEqual(result, frozen)
        self.assertEqual(semantics["rng_relation"], "independent_parent_and_child_streams")
        self.assertEqual(list(semantics["topology_anchors"]), ["thermal_even", "matching_function", "P4_S", "P4_D"])


if __name__ == "__main__":
    unittest.main()
