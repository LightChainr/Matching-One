import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_pivotal_h4_bridge import bridge_metrics  # noqa: E402


class PivotalH4BridgeTests(unittest.TestCase):
    def test_bridge_normalizations(self):
        row = {
            "N": 16,
            "Mbar_prime": 8.0,
            "P4_S_prime": 2.0,
            "P4_D_prime": -1.0,
            "P4_D": 0.5,
        }
        result = bridge_metrics(row)
        self.assertEqual(result["normalized_pivotal_H4"], 0.25)
        self.assertEqual(result["pivotal_H4_scaled"], 0.25 * 16 ** (13 / 8))
        self.assertEqual(result["even_pivotal_H4_scaled"], -2.0)
        self.assertEqual(result["coefficient_ratio"], 0.5)
        self.assertEqual(result["thermal_mass_scaled"], 8.0 / 16 ** (3 / 8))


if __name__ == "__main__":
    unittest.main()
