from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from scripts.primitive_homology_c3_orbits import (
    ROTATION,
    apply,
    build_contract,
    c3_orbits,
    character_inner,
    hexagonal_norm,
    primitive_lines,
    projector_exponents,
    rotate_line,
    spin_charge,
    validate_contract,
)


class PrimitiveHomologyC3OrbitTests(unittest.TestCase):
    def test_checked_in_contract_closes_exactly(self) -> None:
        self.assertEqual(validate_contract(), build_contract())

    def test_rotation_has_order_three_modulo_line_sign(self) -> None:
        for line in primitive_lines(31):
            self.assertEqual(hexagonal_norm(rotate_line(line)), hexagonal_norm(line))
            self.assertEqual(rotate_line(rotate_line(rotate_line(line))), line)
        vector = (2, 1)
        self.assertEqual(apply(ROTATION, apply(ROTATION, apply(ROTATION, vector))), (-2, -1))

    def test_orbits_partition_every_bounded_primitive_line(self) -> None:
        lines = primitive_lines(31)
        flattened = [line for orbit in c3_orbits(31) for line in orbit]
        self.assertEqual(sorted(flattened), lines)
        self.assertEqual(len(flattened), len(set(flattened)))

    def test_character_projectors_are_exactly_orthogonal(self) -> None:
        self.assertEqual(
            [[character_inner(a, b) for b in range(3)] for a in range(3)],
            [[3, 0, 0], [0, 3, 0], [0, 0, 3]],
        )
        self.assertEqual(spin_charge(4), 2)
        self.assertEqual(spin_charge(8), 1)
        self.assertEqual(spin_charge(12), 0)
        self.assertEqual(projector_exponents(2), [0, 1, 2])

    def test_contract_drift_fails_closed(self) -> None:
        frozen = build_contract()
        frozen["rotation_order_on_unoriented_lines"] = 6
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "drift.json"
            path.write_text(json.dumps(frozen), encoding="utf-8")
            with self.assertRaisesRegex(AssertionError, "contract drifted"):
                validate_contract(path)


if __name__ == "__main__":
    unittest.main()
