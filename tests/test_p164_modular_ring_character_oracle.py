from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SCRIPT = ROOT / "scripts" / "p164_modular_ring_character_oracle.py"
SPEC = importlib.util.spec_from_file_location("p164_ring", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class ModularRingCharacterOracleTests(unittest.TestCase):
    def test_low_weight_basis_and_first_degeneracy(self) -> None:
        self.assertEqual(MODULE.modular_basis(2), [])
        self.assertEqual(MODULE.modular_basis(4), [(1, 0)])
        self.assertEqual(MODULE.modular_basis(6), [(0, 1)])
        self.assertEqual(MODULE.modular_basis(8), [(2, 0)])
        self.assertEqual(MODULE.modular_basis(10), [(1, 1)])
        self.assertEqual(MODULE.modular_basis(12), [(0, 2), (3, 0)])

    def test_hex_child_character_is_E4_exponent_mod_three(self) -> None:
        for weight in range(4, 25, 2):
            for a, b in MODULE.modular_basis(weight):
                self.assertEqual((4 * a + 6 * b) % 6, (4 * (a % 3)) % 6)

    def test_numeric_three_child_DFT_is_pure(self) -> None:
        oracle = MODULE.numerical_oracle(80)
        for value in oracle["errors"].values():
            self.assertLess(float(value), 1e-68)
        expected_support = {"E4": 1, "E6": 0, "E4^2": 2}
        for name, support in expected_support.items():
            powers = [float(value) for value in oracle["vectors"][name]["normalized_DFT_power"]]
            self.assertGreater(powers[support], 1 - 1e-15)
            self.assertLess(sum(value for i, value in enumerate(powers) if i != support), 1e-15)

    def test_artifact_is_reproducible(self) -> None:
        expected = json.loads((ROOT / "predictions" / "p164_modular_ring_character_syzygy_20260829.json").read_text())
        self.assertEqual(MODULE.analyze(90), expected)


if __name__ == "__main__":
    unittest.main()
