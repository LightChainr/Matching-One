
from __future__ import annotations
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_axis_pair_annihilator_stable import fit_f_shape, fit_root_power  # noqa: E402


class AxisPairAnnihilatorScorerTests(unittest.TestCase):
    @staticmethod
    def synthetic_rows(q: float = 3.0, w: float = 7.0):
        correction = 0.75
        thermal = 2.0e-13
        pc = 0.59274605079
        root_amplitude = -0.02
        rows = []
        for L in range(10, 19):
            d_q = L ** (-q) - (L - 1) ** (-q)
            d_4 = L ** 4 - (L - 1) ** 4
            rows.append(
                {
                    "pair_L": L,
                    "F_p_ref": correction * d_q + thermal * d_4,
                    "F_coupled_se": 2.0e-8,
                    "annihilator_root": pc + root_amplitude * L ** (-w),
                    "se": {"annihilator_root": 1.0e-10},
                }
            )
        return rows

    def test_true_q3_has_zero_heldout_residual_in_exact_synthetic_data(self) -> None:
        rows = self.synthetic_rows()
        score = fit_f_shape(rows, 3.0, train_max_L=14)
        self.assertLess(score["heldout_chi_square"], 1.0e-16)
        self.assertAlmostEqual(score["parameters"]["correction_amplitude"], 0.75, places=9)
        self.assertAlmostEqual(score["parameters"]["thermal_mistuning_nuisance"], 2.0e-13, places=20)

    def test_wrong_q_is_discriminated_when_noise_is_small(self) -> None:
        rows = self.synthetic_rows()
        right = fit_f_shape(rows, 3.0, train_max_L=14)
        wrong = fit_f_shape(rows, 2.0, train_max_L=14)
        self.assertLess(right["heldout_chi_square"], wrong["heldout_chi_square"])
        self.assertGreater(wrong["heldout_chi_square"], 1.0)

    def test_true_w7_root_crosscheck(self) -> None:
        rows = self.synthetic_rows()
        right = fit_root_power(rows, 7.0, train_max_L=14)
        wrong = fit_root_power(rows, 10.0, train_max_L=14)
        # Exact synthetic data still traverse a double-precision weighted
        # normal-equation solve.  The relevant contract is numerical zero and
        # clear separation from the wrong exponent, not sub-1e-12 chi-square.
        self.assertLess(right["heldout_chi_square"], 1.0e-9)
        self.assertLess(right["heldout_chi_square"], wrong["heldout_chi_square"])
        self.assertAlmostEqual(right["parameters"]["pc"], 0.59274605079, places=13)


if __name__ == "__main__":
    unittest.main()
