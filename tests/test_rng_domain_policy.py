from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from rng_domain_policy import derive_size_seed, domain_record  # noqa: E402


class RngDomainPolicyTests(unittest.TestCase):
    def test_locked_derivation_vector(self):
        self.assertEqual(
            derive_size_seed(2026105001, "P50-fullcurve", 290),
            624074519655001464,
        )

    def test_distinct_sizes_are_domain_separated(self):
        seed65 = derive_size_seed(123, "campaign", 65)
        seed85 = derive_size_seed(123, "campaign", 85)
        self.assertNotEqual(seed65, seed85)
        self.assertEqual(seed65, derive_size_seed(123, "campaign", 65))

    def test_same_size_pair_uses_one_common_field(self):
        record = domain_record(123, "campaign", 65)
        self.assertEqual(record["within_same_size_orientation_pair"], "shared_common_field")
        self.assertEqual(record["cross_size_covariance_contract"], "independent_across_sizes")

    def test_coupling_requires_named_residual_and_full_covariance(self):
        with self.assertRaisesRegex(ValueError, "prespecified residual"):
            domain_record(123, "campaign", 65, mode="intentional_cross_size_coupling")
        record = domain_record(
            123, "campaign", 65, mode="intentional_cross_size_coupling",
            coupled_residual="DeltaM_child-r*DeltaM_parent",
        )
        self.assertEqual(record["effective_seed"], 123)
        self.assertIn("full_cross_size_covariance_required", record["cross_size_covariance_contract"])


if __name__ == "__main__":
    unittest.main()
