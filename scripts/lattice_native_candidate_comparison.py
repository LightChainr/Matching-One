#!/usr/bin/env python3
"""Exact comparison of the four frozen lattice-native candidate values."""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.exact_polynomial_root_certificate import isolate_roots, sturm_sequence, open_root_count
except ModuleNotFoundError:
    from exact_polynomial_root_certificate import isolate_roots, sturm_sequence, open_root_count

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "pslq_search_contract.json"
DEFAULT_OUTPUT = ROOT / "results" / "pslq-lattice-native-candidates" / "latest.json"
SCHEMA = "matching-one/lattice-native-candidate-comparison/v1"
POLYNOMIALS = {
    "kagome-site": [1, 0, -3, 1],
    "three-twelve-site": [1, 0, 0, 0, -3, 0, 1],
    "martini-descendant-golden": [-1, 1, 1],
    "martini-descendant-root2": [-1, 0, 2],
}
WINDOWS = {
    "kagome-site": (Fraction(3, 5), Fraction(7, 10)),
    "three-twelve-site": (Fraction(4, 5), Fraction(5, 6)),
    "martini-descendant-golden": (Fraction(3, 5), Fraction(2, 3)),
    "martini-descendant-root2": (Fraction(7, 10), Fraction(3, 4)),
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


@lru_cache(maxsize=2)
def build_result(contract_path: Path = DEFAULT_CONTRACT) -> dict[str, Any]:
    contract_bytes = contract_path.read_bytes()
    contract = json.loads(contract_bytes)
    library = [row["id"] for row in contract["search_stages"]["lattice_native_candidates"]["library"]]
    _require(library == list(POLYNOMIALS), "lattice-native library/order drift")
    candidates = []
    for candidate_id in library:
        polynomial = [Fraction(value) for value in POLYNOMIALS[candidate_id]]
        lo_window, hi_window = WINDOWS[candidate_id]
        sequence = sturm_sequence(polynomial)
        _require(open_root_count(sequence, lo_window, hi_window) == 1, "candidate window is not isolating")
        roots = isolate_roots(polynomial, lo_window, hi_window, bits=100)
        _require(len(roots) == 1, "candidate root isolation drift")
        lo, hi = roots[0]
        comparisons = []
        for row in contract["intervals"]:
            method_lo, method_hi = Fraction(row["lower"]), Fraction(row["upper"])
            intersection = max(lo, method_lo) <= min(hi, method_hi)
            separation = method_lo - hi if hi < method_lo else lo - method_hi if method_hi < lo else Fraction(0)
            comparisons.append({"interval_id": row["id"], "intersects": intersection, "separation_lower_bound": _text(separation)})
        candidates.append({
            "candidate_id": candidate_id,
            "minimal_polynomial_coefficients_ascending": POLYNOMIALS[candidate_id],
            "isolating_interval": [_text(lo), _text(hi)],
            "isolation_bits": 100,
            "method_comparisons": comparisons,
            "excluded_by_all_method_intervals": all(not row["intersects"] for row in comparisons),
        })
    provenance = contract["provenance"]
    digest = hashlib.sha256((ROOT / provenance["path"]).read_bytes()).hexdigest()
    _require(digest == provenance["sha256"], "provenance digest drift")
    return {
        "schema": SCHEMA,
        "issue": 1,
        "status": "lattice_native_candidate_comparison_complete",
        "contract_sha256": hashlib.sha256(contract_bytes).hexdigest(),
        "provenance_sha256": digest,
        "candidates": candidates,
        "conclusion": {"all_four_candidates_excluded": all(row["excluded_by_all_method_intervals"] for row in candidates)},
        "claim_boundary": {
            "included": "exact algebraic isolation and interval intersection for the four frozen candidates",
            "excluded": "candidate-library expansion, matching-partner independence, near-hit promotion, closed forms, or transcendence",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any], contract_path: Path = DEFAULT_CONTRACT) -> Mapping[str, Any]:
    expected = build_result(contract_path)
    _require(result == expected, "lattice-native result does not exactly reproduce")
    return {"schema": SCHEMA, "status": "valid", "candidates": len(expected["candidates"])}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        print(json.dumps(validate_result(json.loads(args.validate.read_text()), args.contract), indent=2, sort_keys=True)); return 0
    rendered = json.dumps(build_result(args.contract), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True); args.output.write_text(rendered)
    else: print(rendered, end="")
    return 0


if __name__ == "__main__": raise SystemExit(main())
