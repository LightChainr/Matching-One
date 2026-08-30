from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_z5_projective_leg_bivariate_state import (  # noqa: E402
    DEGREE4_HOLDOUT,
    MIXED_COMMUTING_GATE,
    fit_model,
    model_residual,
    recurrence_roots,
)
from score_z5_projective_leg_pair_transfer import CHANNELS  # noqa: E402
from z5_projective_leg_bivariate_mc import GRID, label  # noqa: E402


def synthetic_values(rank: int = 2):
    xroots = (0.62 + 0.04j, -0.11 - 0.31j)[:rank]
    yroots = (-0.11 - 0.31j, 0.62 + 0.04j)[:rank]
    values = {}
    for channel_index, (hand, charge) in enumerate(CHANNELS):
        amplitudes = (1.0 + 0.1j * channel_index, 0.2 - 0.03j * channel_index)[:rank]
        for a, b in GRID:
            value = sum(
                amplitudes[index] * xroots[index] ** a * yroots[index] ** b
                for index in range(rank)
            )
            prefix = f"{label(a, b)}_r{charge}_{hand}_"
            values[prefix + "re"] = value.real
            values[prefix + "im"] = value.imag
    return values


class Z5ProjectiveLegBivariateStateTests(unittest.TestCase):
    def test_durand_kerner_recovers_quadratic_roots(self) -> None:
        expected = (0.6 + 0.1j, -0.2 - 0.3j)
        coefficients = [sum(expected), -expected[0] * expected[1]]
        observed = recurrence_roots(coefficients)
        self.assertLess(min(abs(observed[0] - root) for root in expected), 1e-12)
        self.assertLess(min(abs(observed[1] - root) for root in expected), 1e-12)

    def test_rank_two_axis_fit_predicts_mixed_and_degree_four(self) -> None:
        values = synthetic_values(2)
        model = fit_model(values, 2)
        mixed = model_residual(values, model, MIXED_COMMUTING_GATE)
        heldout = model_residual(values, model, DEGREE4_HOLDOUT)
        self.assertLess(max(abs(value) for value in (*mixed, *heldout)), 1e-10)


if __name__ == "__main__":
    unittest.main()
