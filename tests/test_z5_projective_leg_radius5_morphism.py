from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from z5_projective_leg_radius5_morphism_mc import (  # noqa: E402
    FIELD_ORDER,
    SHELL,
    SOURCE_NEW,
    d4_closure,
    exact_gate,
)


class Z5ProjectiveLegRadius5MorphismTests(unittest.TestCase):
    def test_three_source_endpoints_force_the_full_shell(self) -> None:
        self.assertEqual(d4_closure(SOURCE_NEW), set(SHELL))

    def test_exact_geometry_gate(self) -> None:
        gate = exact_gate()
        self.assertTrue(gate["passed"])
        self.assertEqual(gate["shell_points"], 20)
        self.assertEqual(gate["distinct_parent_vertices"], 20)
        self.assertEqual(gate["complex_rows_per_batch"], 80)
        self.assertEqual(len(FIELD_ORDER), 160)


if __name__ == "__main__":
    unittest.main()
