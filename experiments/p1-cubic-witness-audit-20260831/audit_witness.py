#!/usr/bin/env python3
"""Audit one already selected cubic with exact arithmetic; never run the search."""
from __future__ import annotations

import argparse
import ast
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import reduce
import hashlib
import json
from math import comb, gcd
from pathlib import Path
import platform
import subprocess
import time

COMMIT = "29902bbd4692f789547f252d121858102f1f6650"
COEFFICIENTS = (-63, 40, 65, 79)
ROOT = Path(__file__).resolve().parent
CONTRACT = "analysis/pslq_search_contract.json"
ENGINE = "scripts/degree3_interval_exclusion.py"
IDS = ("jacobsen-2015-eigenvalue", "mertens-2022-p-med", "mertens-2022-p-cell", "yang-zhou-2024-corrected")
PATHS = [CONTRACT, ENGINE] + [f"results/pslq-degree3-{i}/latest.json" for i in IDS]


def text(x):
    x = Fraction(x)
    return f"{x.numerator}/{x.denominator}"


def decimal(x, digits=40):
    x = Fraction(x)
    with localcontext() as ctx:
        ctx.prec = digits
        return str(Decimal(x.numerator) / Decimal(x.denominator))


def rational(x, digits=40):
    return {"exact": text(x), "decimal": decimal(x, digits)}


def poly(x, coefficients=COEFFICIENTS):
    answer = Fraction(0)
    for a in reversed(coefficients):
        answer = answer * x + a
    return answer


def bracket(bits):
    # P is globally increasing, P(0)<0<P(1). All decisions below are exact.
    denominator = 1 << bits
    low, high = 0, denominator
    while high - low > 1:
        mid = (low + high) // 2
        v = poly(Fraction(mid, denominator))
        assert v != 0, "Unexpected dyadic rational root"
        if v < 0:
            low = mid
        else:
            high = mid
    lo, hi = Fraction(low, denominator), Fraction(high, denominator)
    assert poly(lo) < 0 < poly(hi)
    return lo, hi


def divisors(n):
    return [d for d in range(1, abs(n) + 1) if n % d == 0]


def load_sources(repo):
    sources = []
    for path in PATHS:
        target = ROOT / "inputs" / path
        if repo:
            data = subprocess.check_output(["git", "show", f"{COMMIT}:{path}"], cwd=repo)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
        data = target.read_bytes()
        sources.append({"git_commit": COMMIT, "repository_path": path,
                        "snapshot_path": str(target.relative_to(ROOT)), "bytes": len(data),
                        "sha256": hashlib.sha256(data).hexdigest()})
    if not repo:
        previous = json.loads((ROOT / "SOURCES.json").read_text())
        assert previous["files"] == sources, "Input snapshots changed"
    result = {"repository": "https://github.com/LightChainr/Matching-One", "commit": COMMIT,
              "files": sources, "access": "fixed git blobs; no enumeration, network or cloud"}
    (ROOT / "SOURCES.json").write_text(json.dumps(result, indent=2) + "\n")
    return sources


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-repo", type=Path, help="First run only: vendor the six fixed input blobs")
    args = parser.parse_args()
    start = time.monotonic()
    sources = load_sources(args.source_repo)
    contract = json.loads((ROOT / "inputs" / CONTRACT).read_text())
    sha_contract = next(s["sha256"] for s in sources if s["repository_path"] == CONTRACT)
    a0, a1, a2, a3 = COEFFICIENTS
    derivative = (a1, 2 * a2, 3 * a3)
    discriminant = derivative[1] ** 2 - 4 * derivative[2] * derivative[0]
    vertex = Fraction(-derivative[1], 2 * derivative[2])
    minimum_derivative = poly(vertex, derivative)
    assert discriminant < 0 and minimum_derivative > 0
    lo80, hi80 = bracket(80)
    lo600, hi600 = bracket(600)
    assert lo80 < lo600 < hi600 < hi80
    assert decimal(lo600, 150) == decimal(hi600, 150), "Display digits are not certified stable"
    root_mid = (lo600 + hi600) / 2
    candidates = sorted({Fraction(s * n, d) for n in divisors(a0) for d in divisors(a3) for s in (-1, 1)})
    rational_roots = [c for c in candidates if poly(c) == 0]
    assert not rational_roots
    complement = tuple(sum(COEFFICIENTS[k] * comb(k, j) * (-1) ** j
                           for k in range(j, 4)) for j in range(4))
    assert complement == (121, -407, 302, -79)
    # Independent exact substitution check, not a numerical root finder.
    assert all(poly(1 - x) == poly(x, complement) for x in (Fraction(0), Fraction(1, 3), Fraction(1), Fraction(2)))
    rows = []
    for interval in contract["intervals"]:
        source = json.loads((ROOT / "inputs" / f"results/pslq-degree3-{interval['id']}/latest.json").read_text())
        stored = source["interval_result"]
        closest = stored["closest_polynomial"]
        lower, upper = Fraction(interval["lower"]), Fraction(interval["upper"])
        central = Fraction(interval["central_value"])
        pl, pu = poly(lower), poly(upper)
        assert pl < pu < 0 and upper < lo80
        assert source["contract_sha256"] == sha_contract
        assert [stored["lower"], stored["upper"]] == [interval["lower"], interval["upper"]]
        assert closest["coefficients_ascending"] == list(COEFFICIENTS)
        assert [Fraction(x) for x in closest["polynomial_endpoint_values"]] == [pl, pu]
        assert Fraction(closest["minimum_absolute_residual"]) == -pu
        certificate = closest["independent_sturm_certificate"]
        assert [[Fraction(x) for x in pair] for pair in certificate["unit_interval_root_brackets"]] == [[lo80, hi80]]
        assert Fraction(certificate["nearest_root_separation_lower_bound"]) == lo80 - upper
        assert certificate["nearest_root_side"] == "above"
        assert stored["root_containing_polynomials"] == stored["derivative_stationary_fibers"] == 0
        assert stored["primitive_cubics_checked"] == 749507743
        d = root_mid - upper
        rows.append({"interval_id": interval["id"], "lower": text(lower), "upper": text(upper),
                     "width": rational(upper - lower), "root_side": "above",
                     "root_distance_to_interval_80bit_bounds": [text(lo80 - upper), text(hi80 - upper)],
                     "root_distance_to_interval_decimal": decimal(d),
                     "root_minus_central_decimal": decimal(root_mid - central),
                     "signed_residual_interval": [rational(pl), rational(pu)],
                     "minimum_absolute_residual": rational(-pu),
                     "derivative_interval": [rational(poly(lower, derivative)), rational(poly(upper, derivative))],
                     "width_over_root_distance_decimal": decimal((upper - lower) / d),
                     "source_witness_endpoint_residual_and_80bit_certificate_match": True,
                     "reported_full_search": {k: stored[k] for k in ("primitive_cubics_checked", "coefficient_fibers_checked", "root_containing_polynomials", "derivative_stationary_fibers", "excluded")}})
    hull_lower = min(Fraction(i["lower"]) for i in contract["intervals"])
    hull_upper = max(Fraction(i["upper"]) for i in contract["intervals"])
    engine_text = (ROOT / "inputs" / ENGINE).read_text()
    tree = ast.parse(engine_text)
    imports = sorted({n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) and n.module}
                     | {a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names})
    assert not ({"subprocess", "ctypes", "cffi", "cppyy"} & set(imports))
    result = {"schema": "matching-one/cubic-witness-audit/v1", "commit": COMMIT,
              "status": "existing_witness_independently_audited_no_search_replay",
              "coefficients_ascending": list(COEFFICIENTS), "height": max(map(abs, COEFFICIENTS)),
              "primitive_gcd": reduce(gcd, map(abs, COEFFICIENTS)),
              "monotonicity": {"derivative_coefficients_ascending": list(derivative),
                               "derivative_discriminant": discriminant,
                               "derivative_vertex": rational(vertex), "global_minimum_derivative": rational(minimum_derivative),
                               "conclusion": "P strictly increasing on R; exactly one real root"},
              "root": {"decimal_150_digits": decimal(root_mid, 150),
                       "decimal_derivation": "midpoint of exact 600-bit dyadic bracket, rounded to 150 significant digits",
                       "bracket_80bit": [text(lo80), text(hi80)], "bracket_80bit_width": text(hi80 - lo80),
                       "bracket_80bit_endpoint_polynomials": [text(poly(lo80)), text(poly(hi80))],
                       "bracket_600bit": [text(lo600), text(hi600)],
                       "unrounded_midpoint_absolute_error_bound": text((hi600 - lo600) / 2)},
              "rational_root_test": {"candidate_count": len(candidates),
                                     "candidates_and_exact_values": [[text(x), text(poly(x))] for x in candidates],
                                     "rational_roots": [], "irreducible_over_Q": True,
                                     "reason": "primitive cubic with no rational root"},
              "complement": {"coefficients_of_P_1_minus_p_ascending": list(complement),
                             "positive_leading_coefficients_ascending": [-x for x in complement],
                             "height": max(map(abs, complement)), "primitive_gcd": reduce(gcd, map(abs, complement)),
                             "inside_original_height_100_box": False,
                             "root_decimal_100_digits": decimal(1 - root_mid, 100),
                             "claim": "algebraic substitution of the same witness; no independent evidence"},
              "intervals": rows,
              "four_interval_geometry": {"hull_lower": text(hull_lower), "hull_upper": text(hull_upper),
                                         "hull_width": rational(hull_upper - hull_lower),
                                         "hull_width_over_nearest_root_distance_decimal": decimal((hull_upper - hull_lower) / (root_mid - hull_upper)),
                                         "witness_has_no_root_in_hull": True,
                                         "not_a_combined_statistical_interval": True,
                                         "ranking_stability_proved": False,
                                         "reason": "No runner-up gap supplied; witness geometry alone cannot certify winner ranking across the hull"},
              "engine_audit": {"implementation": "Python", "arithmetic": "arbitrary precision integers and fractions.Fraction",
                               "imports": imports, "cpp17_search_engine": False,
                               "objective": "minimum over candidate P of minimum over p in the method interval of abs(P(p))",
                               "tie_break": "lexicographic ascending coefficient tuple",
                               "root_isolation": "only the residual-minimizing selected polynomial is passed to the Sturm root certificate",
                               "globally_smallest_root_distance_certified": False,
                               "execution_receipt_or_CI_environment_in_supplied_inputs": False},
              "boundaries": ["No 749,507,743-candidate enumeration rerun; cached full-search conclusions are reported, not independently regenerated.",
                             "Four method intervals remain separate and quoted uncertainties are not confidence-homogenized.",
                             "No p-values, coefficient physics, closed-form threshold identification, transcendence claim, runner-up gap or globally nearest-root claim.",
                             "Finite degree-3 height-100 exclusions do not settle other degrees/heights or the parent scientific question."]}
    (ROOT / "latest.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    receipt = {"python": platform.python_version(), "executable": __import__("sys").executable,
               "elapsed_seconds": time.monotonic() - start, "new_samples": 0, "enumeration_calls": 0,
               "cloud_calls": 0, "fixed_commit": COMMIT, "audit_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    with (ROOT / "run_receipts.jsonl").open("a") as out:
        out.write(json.dumps(receipt) + "\n")
    print(json.dumps({"root": result["root"]["decimal_150_digits"], "intervals": [{"id": r["interval_id"], "distance": r["root_distance_to_interval_decimal"], "residual": r["minimum_absolute_residual"]["decimal"]} for r in rows], "elapsed_seconds": receipt["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    run()
