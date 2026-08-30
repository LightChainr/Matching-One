from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_p154_phase_e_mixed_plane_pilot as score


def test_mixed_row_algebra() -> None:
    row = {
        "n": "10", "samples": "10", "sum_i0": "2", "sum_i2": "3",
        "sum_k1": "20", "sum_k2": "50",
        "sum_black_axis_pairs": "40", "sum_white_matching_axis_pairs": "60",
        "sum_even_numerator_squared": "1200",
        "sum_i0_even_numerator": "18", "sum_i2_even_numerator": "35",
    }
    values = score.row_values(row)
    assert set(values) == set(score.FIELDS)
    assert values["B"] == 0.25
    assert abs(values["J_top"] - ((10 / 9) * (35 / 400 - 0.3 * 0.25)
                                  - (10 / 9) * (18 / 400 - 0.2 * 0.25))) < 1e-15


def test_protocol_has_all_mixed_rows() -> None:
    manifest = (ROOT / "experiments" / "p154_phase_e_mixed_plane_pilot_20260830.json").read_text()
    for token in ("B^2", "I0*B", "I2*B", "J_top", "J_bulk"):
        assert token in manifest
