from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


def cos4(a: int, b: int) -> Fraction:
    n = a * a + b * b
    return Fraction(a**4 - 6 * a * a * b * b + b**4, n * n)


class MeanJD4Q2ChainPreregTests(unittest.TestCase):
    def test_exact_alternating_geometry_and_frozen_stream(self) -> None:
        manifest = json.loads(
            (ROOT / "analysis" / "mean_jd4_q2_chain_prereg.json").read_text(encoding="utf-8")
        )
        self.assertTrue(manifest["frozen_before_N260_outcome"])
        self.assertEqual(manifest["samples_per_size"], 5_000_000)
        deltas = []
        for row in manifest["sizes"]:
            first, second = row["first"], row["second"]
            delta = cos4(*first) - cos4(*second)
            self.assertEqual(delta, Fraction(row["delta_cos4"]))
            deltas.append(delta)
        self.assertEqual(deltas, [Fraction(1152, 845), Fraction(-1152, 845), Fraction(1152, 845)])
        self.assertEqual(manifest["batches"], 100)
        self.assertEqual(manifest["replica_offset"], 9_300_000_000)


if __name__ == "__main__":
    unittest.main()
