#!/usr/bin/env python3
"""Frozen, exact N25/m64 two-law comparison; no enumeration or sampling.

Run only after committing this package:
    python3 score.py --freeze-commit FULL_SHA --output-dir NEW_DIRECTORY

All numerical decisions use rational bounds. Floating values are display only.
The vendored Interval implementation is loaded only after its provenance gate.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import re
import subprocess
import sys
import time
import types


PACKAGE = Path(__file__).resolve().parent
SOURCE_COMMIT = "a70eeff09f51ce2fa0fea5ae637e9191efbf2e1f"
N, M, BISECTIONS = 25, 64, 160
DELTA = F(1152, 625)
TARGET_PROBABILITY = F(19, 20)
SAMPLING_GATE = 10**9
ROLES = {"axis_histogram", "tilted_histogram", "interval_backend"}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def git(repo, *args):
    return subprocess.run(
        ["git", "-C", str(repo), *args], check=True, capture_output=True
    ).stdout


def provenance_gate(freeze_commit):
    """Read and pin every consumed file before any histogram evaluation."""
    if not re.fullmatch(r"[0-9a-f]{40}", freeze_commit):
        raise ValueError("--freeze-commit must be a full lowercase 40-digit SHA")
    repo = Path(git(PACKAGE, "rev-parse", "--show-toplevel").decode().strip())
    resolved = git(repo, "rev-parse", freeze_commit + "^{commit}").decode().strip()
    if resolved != freeze_commit:
        raise ValueError("freeze argument did not resolve to that exact commit")
    relative_package = PACKAGE.relative_to(repo)
    files = {}
    for name in ("CONTRACT.md", "score.py", "inputs/SOURCES.json"):
        files[name] = (PACKAGE / name).read_bytes()
    sources = json.loads(files["inputs/SOURCES.json"])
    entries = sources["entries"]
    if len(entries) != 3 or {e["role"] for e in entries} != ROLES:
        raise ValueError("exactly the two histograms and one Interval backend are required")
    by_role = {}
    for entry in entries:
        name = entry["path"]
        path = Path(name)
        if path.is_absolute() or ".." in path.parts or path.as_posix() != name:
            raise ValueError("input path must be canonical and package-relative")
        if name in files:
            raise ValueError("duplicate or reserved input path")
        if not (PACKAGE / path).resolve().is_relative_to(PACKAGE):
            raise ValueError("input symlink escapes package")
        if entry["source_commit"] != SOURCE_COMMIT:
            raise ValueError("all inputs must come from the fixed a70eeff0 commit")
        blob = (PACKAGE / path).read_bytes()
        if sha(blob) != entry["sha256"]:
            raise ValueError("input SHA256 mismatch: " + name)
        original = git(repo, "show", SOURCE_COMMIT + ":" + entry["source_path"])
        if original != blob:
            raise ValueError("input differs from declared original Git blob: " + name)
        files[name] = blob
        by_role[entry["role"]] = entry
    frozen_files = {}
    for name, blob in files.items():
        git_path = (relative_package / name).as_posix()
        frozen_blob = git(repo, "show", freeze_commit + ":" + git_path)
        if frozen_blob != blob:
            raise ValueError("package differs from freeze commit: " + git_path)
        frozen_files[name] = {
            "sha256": sha(blob), "bytes": len(blob), "freeze_git_path": git_path,
        }
    receipt = {
        "freeze_commit": freeze_commit,
        "source_commit": SOURCE_COMMIT,
        "package_git_path": relative_package.as_posix(),
        "verified_files": frozen_files,
        "original_sources": entries,
        "verification": "byte equality to freeze and original input Git blobs; SHA256 equality",
    }
    return files, by_role, receipt


def load_backend(data, path):
    # Execute only the already verified byte snapshot, never reread a moving file.
    name = "p337_m64_frozen_interval_backend"
    module = types.ModuleType(name)
    module.__file__ = str(PACKAGE / path)
    sys.modules[name] = module
    exec(compile(data, module.__file__, "exec"), module.__dict__)
    return module


def read_histogram(data):
    reader = csv.DictReader(io.StringIO(data.decode("utf-8")))
    if reader.fieldnames != ["k", "g", "q", "count"]:
        raise ValueError("histogram columns must be k,g,q,count in that order")
    rows, seen = [], set()
    totals = [0] * (N + 1)
    for row in reader:
        k, g, q, count = (int(row[key]) for key in reader.fieldnames)
        if not (0 <= k <= N and g >= 0 and q in (-1, 0, 1) and count > 0):
            raise ValueError("invalid histogram entry")
        if (k, g, q) in seen:
            raise ValueError("duplicate histogram bin")
        seen.add((k, g, q))
        totals[k] += count
        rows.append((k, g, q, count))
    expected = [math.comb(N, k) for k in range(N + 1)]
    if totals != expected or sum(totals) != 2**N:
        raise ValueError("histogram does not contain all binomial(25,K) configurations")
    return rows, {"nonzero_bins": len(rows), "counts_by_K": totals, "total": sum(totals)}


def coefficients(rows, law):
    if law not in ("star", "drop"):
        raise ValueError("only the two frozen laws are permitted")
    c = {key: [F(0)] * (N + 1) for key in ("z", "q", "e", "rank1")}
    for k, g, q, count in rows:
        r = q + 1
        weight = F(count) * F(M)**(-g + (r if law == "drop" else 0))
        c["z"][k] += weight
        c["q"][k] += q * weight
        c["e"][k] += q*q * weight
        c["rank1"][k] += int(q == 0) * weight
    if any(c["e"][k] + c["rank1"][k] != c["z"][k] for k in range(N + 1)):
        raise ArithmeticError("E=1-P(rank1) coefficient identity failed")
    return c


def poly(c, h):
    value = 0
    for coefficient in reversed(c):
        value = value*h + coefficient
    return value


def derivative(c):
    return [k*c[k] for k in range(1, len(c))] or [F(0)]


def convolution(a, b):
    result = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return result


def add(a, b, sign=1):
    result = [F(0)] * max(len(a), len(b))
    for i, x in enumerate(a):
        result[i] += x
    for i, x in enumerate(b):
        result[i] += sign*x
    return result


def primitive_integer_polynomial(c):
    c = list(c)
    while c and c[-1] == 0:
        c.pop()
    if not c:
        raise ArithmeticError("root numerator vanishes identically")
    denominator = math.lcm(*(x.denominator for x in c))
    integers = [int(x*denominator) for x in c]
    divisor = math.gcd(*integers)
    return [x // divisor for x in integers]


def isolate_root(axis, tilted, Interval):
    # 2 Q = qa/Za + qt/Zt; both partition polynomials are positive for h>0.
    numerator = add(convolution(axis["q"], tilted["z"]),
                    convolution(tilted["q"], axis["z"]))
    integer_numerator = primitive_integer_polynomial(numerator)
    signs = [1 if x > 0 else -1 for x in integer_numerator if x]
    variations = sum(a != b for a, b in zip(signs, signs[1:]))
    if variations != 1:
        raise ArithmeticError(
            "STOP: Descartes variation is %d, not one; do not bypass this gate" % variations
        )
    if integer_numerator[0] >= 0 or integer_numerator[-1] <= 0:
        raise ArithmeticError("STOP: expected negative constant and positive leading coefficient")
    lo, hi = F(0), F(1)
    doublings = 0
    while poly(integer_numerator, hi) <= 0:
        # This is endpoint bracketing in h, not a scan of the fixed multiplier m.
        hi *= 2
        doublings += 1
        if doublings > 256:
            raise ArithmeticError("STOP: deterministic h endpoint bracket failed")
    initial_hi = hi
    if not poly(integer_numerator, lo) < 0 < poly(integer_numerator, hi):
        raise ArithmeticError("STOP: initial rational sign bracket failed")
    exact_hit = False
    completed = 0
    for _ in range(BISECTIONS):
        midpoint = (lo + hi) / 2
        value = poly(integer_numerator, midpoint)
        completed += 1
        if value == 0:
            lo = hi = midpoint
            exact_hit = True
            break
        if value < 0:
            lo = midpoint
        else:
            hi = midpoint
    final_signs = (0, 0) if exact_hit else (-1, 1)
    if exact_hit:
        if poly(integer_numerator, lo) != 0:
            raise ArithmeticError("exact rational root check failed")
    elif not poly(integer_numerator, lo) < 0 < poly(integer_numerator, hi):
        raise ArithmeticError("final rational sign bracket failed")
    if lo <= 0:
        raise ArithmeticError("STOP: 160 bisections did not isolate a strictly positive root")
    certificate = {
        "polynomial": "q_axis*Z_tilt + q_tilt*Z_axis",
        "primitive_integer_coefficients_ascending": [str(x) for x in integer_numerator],
        "degree": len(integer_numerator)-1,
        "nonzero_coefficient_signs_ascending": signs,
        "descartes_sign_variations": variations,
        "unique_positive_root": True,
        "partition_denominators_strictly_positive_for_h_positive": True,
        "initial_lower_fraction": "0", "initial_upper_fraction": str(initial_hi),
        "bracket_doublings_from_one": doublings,
        "max_bisections": BISECTIONS, "bisections_completed": completed,
        "exact_rational_root_hit": exact_hit,
        "h_root_lower_fraction": str(lo), "h_root_upper_fraction": str(hi),
        "root_numerator_endpoint_signs": list(final_signs),
        "h_root_midpoint_approx": float((lo + hi) / 2),
    }
    return Interval(lo, hi), certificate


def normalized_observers(c, h):
    z = poly(c["z"], h)
    if z.lo <= 0:
        raise ArithmeticError("partition enclosure must be strictly positive")
    result = {"Z": z}
    z_prime = derivative(c["z"])
    for key in ("q", "e", "rank1"):
        result[key] = poly(c[key], h) / z
        # Cancel algebraically before interval evaluation. This is the exact
        # quotient derivative, including the full geometry-specific normalizer.
        numerator = add(convolution(derivative(c[key]), c["z"]),
                        convolution(c[key], z_prime), sign=-1)
        result[key + "_h"] = poly(numerator, h) / (z*z)
    # Exact coefficient identity independently checks the thermal channel mapping.
    derivative_e = add(convolution(derivative(c["e"]), c["z"]),
                       convolution(c["e"], z_prime), sign=-1)
    derivative_p1 = add(convolution(derivative(c["rank1"]), c["z"]),
                        convolution(c["rank1"], z_prime), sign=-1)
    if any(add(derivative_e, derivative_p1)):
        raise ArithmeticError("E_h=-P1_h polynomial identity failed")
    return result


def bound_sign(serialized):
    lo, hi = F(serialized["lower_fraction"]), F(serialized["upper_fraction"])
    return "positive" if lo > 0 else "negative" if hi < 0 else "unresolved"


def rank1_sample_bound(serialized):
    upper = F(serialized["upper_fraction"])
    if upper < 0:
        raise ArithmeticError("negative rank-one probability upper bound")
    if upper == 0:
        return {"necessary_draws_lower_bound": None, "finite_budget_possible": False,
                "exceeds_one_billion": True,
                "interpretation": "rank-one event has zero probability"}
    bound = math.ceil(TARGET_PROBABILITY / upper)
    return {
        "necessary_draws_lower_bound": bound,
        "finite_budget_possible": True,
        "target_probability_fraction": str(TARGET_PROBABILITY),
        "P1_upper_fraction_used": str(upper),
        "formula": "ceil((19/20)/P1_upper)",
        "exceeds_one_billion": bound > SAMPLING_GATE,
        "interpretation": "necessary only for seeing at least one rank-one draw; not sufficient for U",
    }


def evaluate_law(law, rows, backend):
    coeffs = {g: coefficients(rows[g], law) for g in ("axis", "tilted")}
    root, certificate = isolate_root(coeffs["axis"], coeffs["tilted"], backend.Interval)
    obs = {g: normalized_observers(coeffs[g], root) for g in ("axis", "tilted")}
    slope = (obs["axis"]["q_h"] + obs["tilted"]["q_h"]) / 2
    slope_json = backend.interval_json(slope)
    if slope.lo <= 0 or bound_sign(slope_json) != "positive":
        raise ArithmeticError("STOP: fixed-precision rational slope enclosure is not strictly positive")
    numerator = (obs["axis"]["e_h"] - obs["tilted"]["e_h"]) / DELTA
    observer = numerator / slope
    observer_json = backend.interval_json(observer)
    cells = {}
    for geometry in ("axis", "tilted"):
        item = obs[geometry]
        probability = item["rank1"]
        if probability.lo < 0 or probability.hi > 1:
            raise ArithmeticError("STOP: fixed-precision P1 enclosure is outside [0,1]")
        p1 = backend.interval_json(probability)
        cells[geometry] = {
            "q": backend.interval_json(item["q"]),
            "E": backend.interval_json(item["e"]),
            "q_h": backend.interval_json(item["q_h"]),
            "E_h": backend.interval_json(item["e_h"]),
            "P_rank1": p1,
            "ordinary_sampling_necessary_bound": rank1_sample_bound(p1),
        }
    return {
        "law": law, "rank_coefficient": 0 if law == "star" else 1,
        "weights": "count*h^K*m^(-g)" if law == "star" else "count*h^K*m^(-g+q+1)",
        "root_certificate": certificate,
        "pooled_Q_at_root": backend.interval_json((obs["axis"]["q"] + obs["tilted"]["q"]) / 2),
        "positive_slope_D": slope_json,
        "thermal_difference_over_Delta": backend.interval_json(numerator),
        "U_over_A25": observer_json,
        "U_sign": bound_sign(observer_json),
        "geometry": cells,
        "normalization": "separate full Z for each geometry at this law's own pooled root",
        "exact_observer_identities": ["E=1-P_rank1", "E_h=-P_rank1_h"],
    }


def render_report(result):
    lines = [
        "# Fixed N25, m64 two-law result", "",
        "Freeze: `" + result["provenance"]["freeze_commit"] + "`.", "",
        "Primary decision: **" + result["primary_decision"] + "**.", "",
        "The displayed bounds are outward-rounded rational enclosures; numerical midpoints do not decide signs. "
        "A25=25^(13/8)/2 is positive, so U and U/A25 have the same sign. Each law uses its own pooled root.", "",
        "| Law | U/A25 lower | U/A25 upper | Sign |", "| --- | --- | --- | --- |",
    ]
    for law in ("star", "drop"):
        row = result["laws"][law]
        u = row["U_over_A25"]
        lines.append(f"| {law} | {u['lower_fraction']} | {u['upper_fraction']} | {row['U_sign']} |")
    lines += ["", "| Law / geometry | P(rank1) lower | P(rank1) upper | Necessary draws lower bound |",
              "| --- | --- | --- | --- |"]
    for law in ("star", "drop"):
        for geometry in ("axis", "tilted"):
            row = result["laws"][law]["geometry"][geometry]
            p1 = row["P_rank1"]
            bound = row["ordinary_sampling_necessary_bound"]["necessary_draws_lower_bound"]
            lines.append(f"| {law} / {geometry} | {p1['lower_fraction']} | {p1['upper_fraction']} | {bound} |")
    lines += ["", "Resource gate: **" + result["resource_decision"] + "**.", "",
              "The union bound uses ceil((19/20)/P1_upper). It is necessary for a 95% chance of seeing even one "
              "rank-one draw, regardless of dependence. It does not give a sufficient budget to estimate U or "
              "a lower bound for importance, conditional or other variance-reduced estimators. No wall-clock estimate follows.", "",
              "The JSON includes both primitive integer root numerators, Descartes certificates, exact root brackets, "
              "positive slope enclosures, full-normalization observer derivatives, and all four population cells. "
              "No enumeration, simulation, cloud job, source fit or other coupling point was evaluated.", "",
              "This calculation is a deterministic consequence of the existing exact histogram. It is not independent evidence, "
              "a continuum or large-N result, a uniform remainder bound, or a homogeneous continuation theorem. "
              "Failure of finite-point separation would not refute the proved eventual asymptotic signs.", "",
              "Reproduce from the frozen checkout (standard-library Python >=3.10):", "", "```sh",
              "python3 experiments/p337-finite-law-window-20260831/score.py \\",
              "  --freeze-commit " + result["provenance"]["freeze_commit"] + " \\",
              "  --output-dir /tmp/p337-m64-fresh-output", "```", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze-commit", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.expanduser().resolve()
    if output.exists():
        raise FileExistsError("output directory must not already exist: " + str(output))
    started = time.monotonic()
    started_utc = datetime.now(timezone.utc).isoformat()
    files, roles, provenance = provenance_gate(args.freeze_commit)
    backend_entry = roles["interval_backend"]
    backend = load_backend(files[backend_entry["path"]], backend_entry["path"])
    rows, validation = {}, {}
    for geometry, role in (("axis", "axis_histogram"), ("tilted", "tilted_histogram")):
        rows[geometry], validation[geometry] = read_histogram(files[roles[role]["path"]])
    laws = {law: evaluate_law(law, rows, backend) for law in ("star", "drop")}
    star_sign, drop_sign = laws["star"]["U_sign"], laws["drop"]["U_sign"]
    if star_sign == "negative" and drop_sign == "positive":
        primary = "certified_Ustar_negative_Udrop_positive_at_fixed_m64"
    elif star_sign == "positive" or drop_sign == "negative":
        primary = "inconsistent_with_the_fixed_finite_point_sign_prediction"
    else:
        primary = "fixed_finite_point_sign_comparison_unresolved"
    resource_gate = any(
        laws[law]["geometry"][geometry]["ordinary_sampling_necessary_bound"]["exceeds_one_billion"]
        for law in ("star", "drop") for geometry in ("axis", "tilted")
    )
    result = {
        "schema": "p337_finite_law_window_m64_v1",
        "provenance": provenance,
        "fixed_parameters": {"N": N, "m": M, "DeltaCos4_fraction": str(DELTA),
                             "max_root_bisections": BISECTIONS,
                             "geometry": {"axis": [5, 0], "tilted": [4, 3]}},
        "histogram_validation": validation,
        "laws": laws,
        "primary_decision": primary,
        "resource_decision": ("do_not_promote_ordinary_unconditional_sampling_to_P0_new_estimator_required"
                              if resource_gate else "one_billion_necessary_bound_gate_not_triggered_not_a_feasibility_proof"),
        "ordinary_sampling_launched": False,
        "scope": "one finite point and two fixed laws; exact existing counts; no new population",
    }
    # Prepare every output in memory first; a failed mathematical gate publishes no result.
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    report = render_report(result)
    receipt = {
        "freeze_commit": args.freeze_commit, "started_utc": started_utc,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.monotonic() - started,
        "python": sys.version, "platform": platform.platform(),
        "argv": sys.argv,
        "latest_json_sha256": sha(serialized.encode()),
        "REPORT_md_sha256": sha(report.encode()),
        "input_and_frozen_file_checks": provenance["verified_files"],
        "evaluated_multipliers": [M], "evaluated_laws": ["star", "drop"],
        "enumerations": 0, "random_draws": 0,
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / "latest.json").write_text(serialized, encoding="utf-8")
    (output / "REPORT.md").write_text(report, encoding="utf-8")
    (output / "run.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output_dir": str(output), "primary_decision": primary,
                      "resource_decision": result["resource_decision"]}, indent=2))


if __name__ == "__main__":
    main()
