from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SCRIPT = ROOT / "scripts" / "p161_g2_root_geometry_oracle.py"
SPEC = importlib.util.spec_from_file_location("p161_g2", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class G2RootGeometryOracleTests(unittest.TestCase):
    def test_first_C4_allowed_exact_hex_survivor_is_E6_squared(self) -> None:
        survivors = [row for row in MODULE.c4_allowed_ring(36) if row["survives_exact_hex"]]
        self.assertEqual(survivors[0]["weight"], 12)
        self.assertEqual(survivors[0]["E4_power"], 0)
        self.assertEqual(survivors[0]["E6_power"], 2)
        self.assertEqual(survivors[0]["root_N_power"], 6)

    def test_C4_allowed_hex_survivors_have_weights_multiple_of_twelve(self) -> None:
        for row in MODULE.c4_allowed_ring(60):
            if row["survives_exact_hex"]:
                self.assertEqual(row["weight"] % 12, 0)

    def test_numerical_three_geometry_and_child_oracles(self) -> None:
        errors = MODULE.numerical_oracle(80)["errors"]
        for value in errors.values():
            self.assertLess(float(value), 1e-68)

    def test_artifact_is_reproducible(self) -> None:
        expected = json.loads((ROOT / "predictions" / "p161_g2_root_geometry_no_go_20260829.json").read_text())
        self.assertEqual(MODULE.analyze(90), expected)


if __name__ == "__main__":
    unittest.main()
