import sys
from pathlib import Path
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from p398_width8_projected_memory import build_result


class ProjectedMemory(unittest.TestCase):
    def test_finite_memory_equation_and_crossing(self):
        result=build_result()
        for row in result["ray_rows"]:
            check=row["volterra_at_half"]
            self.assertAlmostEqual(check["derivative_from_identity"],check["derivative_from_spectrum"],places=10)
        self.assertGreater(result["ray_rows"][0]["feedback_vs_force_variance_ratio"],3)
        crossing=result["crossing"]
        self.assertLess(abs(crossing["integrated_decay_difference_at_actual_crossing"]),1e-10)
        self.assertLess(crossing["instantaneous_decay_crossings"][0],crossing["actual_normalized_ray_crossing"])


if __name__=="__main__":
    unittest.main()
