from fractions import Fraction
import json
from pathlib import Path

from scripts.p14_four_terminal_balance_roots import (
    ROOT,
    bernstein_to_power,
    build_result,
    eval_poly,
    isolate_roots,
)


def test_bernstein_conversion_preserves_values():
    counts = [1, 4, 0, 0, 0]
    power = bernstein_to_power(counts)
    for point in (Fraction(0), Fraction(1, 3), Fraction(1)):
        direct = sum(
            count * point**k * (1 - point) ** (4 - k)
            for k, count in enumerate(counts)
        )
        assert eval_poly(power, point) == direct


def test_sturm_isolates_simple_and_exact_roots():
    # (p-1/4)(p-1/2)(p-3/4)
    poly = [Fraction(-3, 32), Fraction(11, 16), Fraction(-3, 2), Fraction(1)]
    roots = isolate_roots(poly)
    assert len(roots) == 3
    assert roots[1][0] <= Fraction(1, 2) <= roots[1][1]
    assert roots[0][0] <= Fraction(1, 4) <= roots[0][1]
    assert roots[2][0] <= Fraction(3, 4) <= roots[2][1]


def test_frozen_corpus_scores_all_candidates():
    source = ROOT / "results" / "terminal-reliability" / "bounded-four-terminal-corpus.json"
    result = build_result(source)
    raw = json.loads(source.read_text(encoding="utf-8"))
    assert result["candidate_count"] == len(raw["candidates"]) == 27
    assert result["candidate_with_open_unit_root_count"] <= 27
    assert all(row["balance_power_coefficients_low_to_high"] for row in result["candidates"])
