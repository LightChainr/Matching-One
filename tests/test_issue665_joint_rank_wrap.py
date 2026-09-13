#!/usr/bin/env python3
"""#665 joint rank x wrap-label census tests (axis L=2,3 — seconds)."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from issue665_joint_rank_wrap_census import census  # noqa: E402


def _tuples(result: dict) -> dict:
    return {
        (r["k"], r["r_black"], r["r_white"], r["wrap_black"], r["wrap_white"],
         r["label5_black"], r["label5_white"]): r["count"]
        for r in result["joint"]
    }


def test_axis_l2_tripwires_and_support() -> None:
    result = census(2)
    assert result["tripwire_bernstein"] == "pass"
    assert result["tripwire_rb_plus_rw_eq_2"] == "pass"
    # rank-pair totals must be exactly the three site classes
    assert result["rank_pair_totals"] == {"0,2": 7, "1,1": 4, "2,0": 5}


def test_axis_l3_both_same_is_not_rank2() -> None:
    result = census(3)
    assert result["tripwire_bernstein"] == "pass"
    assert result["tripwire_rb_plus_rw_eq_2"] == "pass"
    # PR #653 published rank-pair totals, reproduced
    assert result["rank_pair_totals"] == {"0,2": 259, "1,1": 162, "2,0": 91}
    # both-same splits: rank-1 spirals (6) + rank-2 crosses (91) on the black side
    assert result["five_name_by_rank_black"]["both-same"] == [0, 6, 91]
    # the question cell both-same x both-same is PURE rank-1 x rank-1 (spirals)
    assert result["both-same_x_both-same_by_rank_pair"] == {"1,1": 6}
    # coarse both x both coincides with it here (no rank-2 x rank-2 at L=3)
    assert result["coarse-both_x_coarse-both_by_rank_pair"] == {"1,1": 6}
    # both-two has zero mass at every k
    joint = _tuples(result)
    assert not any(k[5] == "both-two" or k[6] == "both-two" for k in joint)
    # D-carrying cells are exactly neither x both (rank 0,2) and both x neither
    # (rank 2,0) — rank-2 sits in the exclusive-cross mismatch, never in bs x bs
    d_cells = {k: v for k, v in joint.items()
               if (k[3] == "neither") != (k[4] == "neither")}
    assert all(k[1:3] == (0, 2) for k in d_cells if k[3] == "neither")
    assert all(k[1:3] == (2, 0) for k in d_cells if k[4] == "neither")


if __name__ == "__main__":
    test_axis_l2_tripwires_and_support()
    test_axis_l3_both_same_is_not_rank2()
    print("ALL_CHECKS_PASS")
