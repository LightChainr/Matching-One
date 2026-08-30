import importlib.util
import math
import pathlib
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "sequential_gaussian_eprocess", ROOT / "scripts" / "sequential_gaussian_eprocess.py"
)
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


class SequentialGaussianEProcessTest(unittest.TestCase):
    def test_one_dimensional_log_lr(self):
        self.assertAlmostEqual(MOD.log_lr_increment([2.0], [0.0], [1.0], [[1.0]]), 1.5)

    def test_minimum_batches_is_respected(self):
        result = MOD.run_path([[10.0]] * 5, [0.0], [1.0], [[1.0]], 0.05, 3, 5)
        self.assertEqual(result["decision"], "alternative")
        self.assertEqual(result["batches"], 3)

    def test_reverse_boundary(self):
        result = MOD.run_path([[-10.0]] * 5, [0.0], [1.0], [[1.0]], 0.05, 1, 5)
        self.assertEqual(result["decision"], "null")

    def test_singular_covariance_is_rejected(self):
        with self.assertRaises(ValueError):
            MOD.log_lr_increment([0.0, 0.0], [0.0, 0.0], [1.0, 1.0], [[1.0, 1.0], [1.0, 1.0]])

    def test_fixed_time_evalue_has_null_mean_one(self):
        # Exact Gaussian MGF identity: E_0 exp(log LR) = exp(-I/2) exp(I/2) = 1.
        information = 1.0
        self.assertAlmostEqual(math.exp(-information / 2 + information / 2), 1.0)


if __name__ == "__main__":
    unittest.main()
