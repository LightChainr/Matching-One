
from __future__ import annotations
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from observer_bandwidth_semantic_gate import build_report, classify  # noqa: E402


class ObserverBandwidthSemanticGateTests(unittest.TestCase):
    def test_valid_product_and_slice_contracts_are_distinct(self) -> None:
        product = classify(
            {
                "base_measure": "bernoulli_product",
                "clock": "product_resample",
                "complete_source_transformed": True,
            }
        )
        slice_result = classify(
            {
                "base_measure": "fixed_k_slice",
                "clock": "occupied_empty_swap",
                "n": 8,
                "k": 4,
                "complete_source_transformed": True,
            }
        )
        self.assertEqual(product["theorem"], "product_walsh")
        self.assertEqual(slice_result["theorem"], "johnson_slice")

    def test_palm_and_adaptive_contracts_fail_closed(self) -> None:
        result = classify(
            {
                "base_measure": "bernoulli_product",
                "clock": "marked_birth",
                "marked_birth": True,
                "adaptive_mark": True,
                "complete_source_transformed": True,
            }
        )
        self.assertFalse(result["accepted"])
        self.assertIn("palm_conditioning_not_removed", result["reasons"])
        self.assertIn("adaptive_mark_degree_not_derived", result["reasons"])

    def test_partial_source_and_degenerate_quotient_fail_closed(self) -> None:
        result = classify(
            {
                "base_measure": "fixed_k_slice",
                "clock": "occupied_empty_swap",
                "n": 8,
                "k": 4,
                "complete_source_transformed": False,
                "degenerate_quotient": True,
            }
        )
        self.assertEqual(
            result["reasons"],
            [
                "complete_source_transform_missing",
                "degenerate_quotient_incidence_missing",
            ],
        )

    def test_slice_endpoints_require_explicit_handling(self) -> None:
        result = classify(
            {
                "base_measure": "fixed_k_slice",
                "clock": "occupied_empty_swap",
                "n": 8,
                "k": 0,
                "complete_source_transformed": True,
            }
        )
        self.assertEqual(result["reasons"], ["slice_endpoint_handling_missing"])

    def test_checked_report_is_exactly_reproducible(self) -> None:
        manifest = json.loads(
            (ROOT / "analysis/observer_bandwidth_semantic_gate_manifest.json").read_text(
                encoding="utf-8"
            )
        )
        checked = json.loads(
            (ROOT / "results/observer-bandwidth-semantic-gate/latest.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(checked, build_report(manifest))
        self.assertEqual(checked["fixture_count"], 7)
        self.assertEqual(
            [row["status"] for row in checked["p267_layers"]],
            ["reject", "accept_johnson_slice", "reject", "reject"],
        )


if __name__ == "__main__":
    unittest.main()
