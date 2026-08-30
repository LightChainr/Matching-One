import json
from pathlib import Path
import unittest

from p334_tm_typed_rewrite_system import build_result


ROOT = Path(__file__).resolve().parents[1]


class P334TMTypedRewriteSystemTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result()

    def test_local_square_and_Alexander_types_are_finite(self):
        theorem = self.result["finite_local_type_theorem"]
        self.assertEqual(theorem["patterns"]["D"]["exit_pattern"], "1222")
        self.assertEqual(
            theorem["patterns"]["D"]["Alexander_dual_pattern"], "0001"
        )
        self.assertEqual(theorem["patterns"]["Y"]["exit_pattern"], "1112")
        self.assertEqual(theorem["unique_negative_product_type"], "D x F")

    def test_canonical_rewrite_closes_every_bounded_row(self):
        audit = self.result["bounded_audit"]
        self.assertEqual(audit["line_layer_rows"], 984)
        self.assertEqual(audit["rewrite_closed_rows"], 984)
        self.assertEqual(audit["rewrite_open_rows"], 0)
        self.assertEqual(audit["unmatched_hard_tokens"], 0)
        self.assertEqual(
            audit["hard_tokens"],
            audit["mixed_tokens_used"] + audit["synergy_tokens_used"],
        )

    def test_cover_cone_rays_have_witnesses_and_rules(self):
        rays = self.result["cover_cone_rays"]
        self.assertEqual(len(rays), 9)
        self.assertEqual(sum(row["on_exact_lower_hull"] for row in rays), 4)
        self.assertEqual(
            sum("then R_Y" in row["canonical_rule"] for row in rays), 3
        )
        self.assertTrue(all(row["minimal_quotient_witness"] for row in rays))

    def test_synergy_rescue_has_four_signatures(self):
        audit = self.result["bounded_audit"]
        self.assertEqual(audit["mixed_only_rows"], 968)
        self.assertEqual(audit["mixed_then_synergy_rows"], 16)
        rescue = self.result["synergy_rescue_classification"]
        self.assertEqual(rescue["signature_count"], 4)
        self.assertEqual(rescue["maximum_synergy_pool_fraction_used"], "133/2880")

    def test_checked_artifact_reproduces(self):
        checked = json.loads(
            (ROOT / "results/p334-tm-typed-rewrite-system/latest.json").read_text()
        )
        self.assertEqual(checked, self.result)


if __name__ == "__main__":
    unittest.main()
