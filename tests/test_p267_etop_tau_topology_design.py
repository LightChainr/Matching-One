from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from check_p267_etop_tau_topology_design import verify  # noqa: E402


class P267EtopTauTopologyDesignTests(unittest.TestCase):
    def test_exact_factorial_and_minimality(self) -> None:
        result = verify()
        self.assertEqual(result["minimal_determinant"], 50)
        self.assertEqual(result["smith"]["cyclic"], [1, 50])
        self.assertEqual(result["delta_chi4"], {"i": "1152/625", "2i": "-1152/625"})


if __name__ == "__main__":
    unittest.main()
