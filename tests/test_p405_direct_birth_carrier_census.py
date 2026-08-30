from fractions import Fraction
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p405_direct_birth_carrier_census import (  # noqa: E402
    build_result,
    carrier_census,
    priority_mass,
)


def test_priority_mass_uses_uniform_permutation_beta_weight():
    counts = [0, 3, 0, 5]
    assert priority_mass(counts, 4) == Fraction(3, 12) + Fraction(5, 4)


def test_n9_carrier_split_closes_total_probability():
    row = carrier_census(3, 0)
    assert (row["direct_edges"], row["theta_edges"], row["figure_eight_edges"]) == (
        45,
        36,
        9,
    )
    assert Fraction(row["theta_probability"]) + Fraction(
        row["figure_eight_probability"]
    ) == Fraction(3, 35)


def test_full_reference_census_and_rows_close():
    result = build_result()
    assert result["decision"] == "theta_and_figure_eight_priority_masses_separated_exactly"
    for geometry in result["geometries"]:
        assert geometry["direct_edges"] == geometry["theta_edges"] + geometry["figure_eight_edges"]
        row_edges = sum(row["direct_edges"] for row in geometry["rows_by_predecessor_size"])
        assert row_edges == geometry["direct_edges"]
