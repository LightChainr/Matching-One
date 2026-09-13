"""Regression lock: the contrast solver must reproduce every published nu control.

The heavy builder (scripts/cylinder_winding_intensity_fast.cpp) is not run here.
This test drives scripts/winding_prefactor_contrast.py against the reward-preserving
lumped tables already committed in
results/geometric-consistency/cylinder-winding-intensity.json, and requires the
recomputed nu to equal the committed exact rational at every control point.
"""
from __future__ import annotations
import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import winding_prefactor_contrast as W  # noqa: E402

REFERENCE = ROOT / "results" / "geometric-consistency" / "cylinder-winding-intensity.json"


class ContrastSolverTests(unittest.TestCase):
    def setUp(self) -> None:
        self.reference = json.loads(REFERENCE.read_text(encoding="utf-8"))

    @staticmethod
    def _aggregate(reduced_table: list) -> list:
        """Reference rows are one entry per mask; aggregate them by popcount."""
        width = len(reduced_table[0]).bit_length() - 1
        out = []
        for row in reduced_table:
            per_m = [dict() for _ in range(width + 1)]
            for mask, (nb, rw) in enumerate(row):
                key = (nb, rw)
                bucket = per_m[bin(mask).count("1")]
                bucket[key] = bucket.get(key, 0) + 1
            out.append([[[nb, rw, c] for (nb, rw), c in sorted(d.items())] for d in per_m])
        return out

    def _as_model(self, model: dict) -> dict:
        """Present a committed model in the shape winding_build emits."""
        return {"width": model["width"], "matching": model["graph"] == "matching",
                "states": model["frontier_states"],
                "rows": self._aggregate(model["reduced_table"])}

    def test_every_published_control_is_reproduced_exactly(self) -> None:
        checked = 0
        for model in self.reference["models"]:
            width = model["width"]
            self.assertEqual(len(model["reduced_table"]), model["reward_lumps"],
                             f"width {width}: lump count disagrees with the header")
            shaped = self._as_model(model)
            for control in model["point_controls"]:
                p = Fraction(control["p"])
                K, g = W.load_counts(shaped, p)
                pi = W.pi_exact_dense(K, g, len(shaped["rows"]))
                nu, _resid, _bound, _gmax, _delta = W.certify(K, g, pi, width, p)
                self.assertEqual(
                    nu, Fraction(control["intensity"]),
                    f"{model['graph']} width {width} p={control['p']}: "
                    f"solver gave {nu}, committed control is {control['intensity']}")
                checked += 1
        self.assertEqual(checked, 18, "expected the full 3x2x3 control grid")

    def test_empty_row_resets_to_the_empty_state(self) -> None:
        """delta = (1-p)^w is the uniform reset probability the certificate divides by."""
        for model in self.reference["models"]:
            width = model["width"]
            rows = model["reduced_table"]
            for control in model["point_controls"]:
                p = Fraction(control["p"])
                shaped = self._as_model(model)
                K, _g = W.load_counts(shaped, p)
                # the mass that lands on lump 0 in the first step from any state is >= (1-p)^w
                reset = (1 - p) ** width
                for i in range(len(rows)):
                    mass = sum(v for j, v in K[i].items() if j == 0)
                    self.assertGreaterEqual(mass, reset - Fraction(1, 10 ** 30))

    def test_nu_is_strictly_between_zero_and_one(self) -> None:
        for model in self.reference["models"]:
            shaped = self._as_model(model)
            for control in model["point_controls"]:
                p = Fraction(control["p"])
                K, g = W.load_counts(shaped, p)
                pi = W.pi_exact_dense(K, g, len(shaped["rows"]))
                nu, _r, _b, _gm, _d = W.certify(K, g, pi, model["width"], p)
                self.assertGreater(nu, 0)
                self.assertLess(nu, 1)


if __name__ == "__main__":
    unittest.main()
