import sys
from pathlib import Path
import unittest

import mpmath as mp


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from p321_aspect_ratio_q4_oracle import (  # noqa: E402
    e4_imaginary,
    oracle_payload,
    parse_positive_rational,
)


class P321AspectRatioQ4OracleTests(unittest.TestCase):
    def test_square_cm_value(self):
        with mp.workdps(100):
            observed = e4_imaginary(mp.mpf(1), dps=80).value
            exact = 3 * mp.gamma(mp.mpf(1) / 4) ** 8 / (64 * mp.pi**6)
            self.assertLess(abs(observed - exact), mp.mpf("1e-75"))

    def test_hecke_ratio_at_two_i(self):
        with mp.workdps(100):
            square = e4_imaginary(mp.mpf(1), dps=80).value
            rectangle = e4_imaginary(mp.mpf(2), dps=80).value
            self.assertLess(abs(rectangle / square - mp.mpf(11) / 16), mp.mpf("1e-75"))

    def test_registered_curve_and_area_conversion(self):
        payload = oracle_payload(dps=70)
        records = payload["records"]
        self.assertEqual([row["rho"] for row in records], ["1", "16/9", "9/4", "4", "9"])
        self.assertEqual(records[0]["width_C_over_square_C"], "1.0")
        for row in records:
            rho = parse_positive_rational(row["rho"])
            width = mp.mpf(row["width_C_over_square_C"])
            area = mp.mpf(row["equal_area_C_over_square_C"])
            self.assertLess(abs(area - rho**2 * width), mp.mpf("1e-60"))
        self.assertLess(abs(mp.mpf(records[-1]["distance_from_cylinder_cusp"])), mp.mpf("7e-23"))

    def test_tail_certificate_dominates_truncation_error(self):
        with mp.workdps(110):
            certified = e4_imaginary(mp.mpf(1), dps=70)
            q = mp.exp(-2 * mp.pi)
            reference = 1 + 240 * mp.fsum(
                sum(d**3 for d in range(1, n + 1) if n % d == 0) * q**n
                for n in range(1, certified.terms + 30)
            )
            self.assertLessEqual(abs(reference - certified.value), certified.absolute_tail_bound)


if __name__ == "__main__":
    unittest.main()
