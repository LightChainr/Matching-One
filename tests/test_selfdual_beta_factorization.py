import importlib.util
import pathlib
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "selfdual_beta_factorization",
    ROOT / "scripts" / "selfdual_beta_factorization.py",
)
MOD = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MOD)


class SelfDualBetaFactorizationTests(unittest.TestCase):
    def test_n10_minimal_degree_beta(self):
        h = [1, 2, 6]
        self.assertEqual(h, MOD.central_binomial_prefix(3))
        self.assertEqual(
            MOD.expand_m_from_h(h),
            [-1, 0, 0, 20, -30, 12],
        )

    def test_n26_symmetric_quotient(self):
        h = [1, 2, 6, 20, 70, 96, 170, 260, 260, 78]
        self.assertEqual(h[:5], MOD.central_binomial_prefix(5))
        self.assertEqual(MOD.minimum_success_count_from_h(5, h[5]), 78)
        self.assertTrue(all(c > 0 for c in h))
        self.assertEqual(
            MOD.expand_m_from_h(h),
            [
                -1, 0, 0, 0, 0, 156, -338, 260, -260, -338,
                1144, 3536, -13702, 15628, -3016, -10088,
                11492, -5798, 1482, -156,
            ],
        )

    def test_beta5_minimum_count_when_degree_is_minimal(self):
        # If H stops at degree s-1, c_s=0 and self-duality/minimum support
        # force the majority-core minimum success count.
        self.assertEqual(MOD.minimum_success_count_from_h(5, 0), 126)


if __name__ == "__main__":
    unittest.main()
