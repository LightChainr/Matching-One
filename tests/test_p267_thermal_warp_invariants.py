from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from p267_thermal_warp_invariants import (certificate, bernstein_to_power,
    clock_moments, density_pullback, deriv, integral, mul, omega, primitive,
    scale, sub, supplied_area_gate)


class WarpInvariantTests(unittest.TestCase):
    def test_exact_counterexamples(self):
        data = certificate()
        self.assertEqual(data["zero_area_but_singular_example"]["flux_at_zero"], "-1/81")
        self.assertEqual(data["finite_TV_obstruction"]["TV_gap"], "17/72")
        self.assertTrue(data["finite_warp_not_linear_generator"]["finite_density_warp_verified"])

    def test_common_warp_area_and_negative_overall_amplitude(self):
        a, e, phi = [F(1), F(1)], [F(1), F(0), F(1)], [F(0), F(5, 4), F(-1, 4)]
        alpha = F(-2, 7)
        ua, ue = scale(density_pullback(a, phi), alpha), scale(density_pullback(e, phi), alpha)
        self.assertEqual(omega(ua, ue), alpha**2*omega(a, e))
        out = supplied_area_gate({"D_A": a, "D_E": e, "U_A": ua, "U_E": ue, "alpha": alpha})
        self.assertEqual(out["common_density_warp_area_null"], "0")

    def test_area_first_variation_with_boundary_mass(self):
        a, e, ra, re = [F(1), F(1)], [F(1), F(0), F(1)], [F(2), F(-1)], [F(3), F(1)]
        ha, he, fa, fe = map(primitive, [ra, re, a, e])
        variation = integral(sub(sub(mul(ha, e), mul(he, a)), sub(mul(fe, ra), mul(fa, re))))
        formula = 2*integral(sub(mul(ha, e), mul(he, a)))+integral(a)*integral(re)-integral(e)*integral(ra)
        self.assertEqual(variation, formula)

    def test_bernstein_input(self):
        self.assertEqual(bernstein_to_power([0, F(1, 2), 0]), [F(0), F(1), F(-1)])

    def test_independent_amplitudes_and_clock_moment_area_equivalence(self):
        a, e, phi = [F(1), F(1)], [F(-1, 3), F(1)], [F(0), F(5, 4), F(-1, 4)]
        ra, re = F(-2, 7), F(3, 5)
        ua, ue = scale(density_pullback(a, phi), ra), scale(density_pullback(e, phi), re)
        out = supplied_area_gate({"D_A": a, "D_E": e, "U_A": ua, "U_E": ue, "r_A": ra, "r_E": re})
        self.assertEqual(out["clock_moment_nulls_m0_to_m6"], ["0"]*7)
        # Perturb E with a zero-area term: identity holds also away from the null.
        ue = sub(ue, [F(-1, 2), F(1)])
        jgap = clock_moments(ua, ue, 1)[1]-re*clock_moments(a, e, 1)[1]
        self.assertEqual(2*ra*integral(a)*jgap, omega(ua, ue)-ra*re*omega(a, e))


if __name__ == "__main__":
    unittest.main()
