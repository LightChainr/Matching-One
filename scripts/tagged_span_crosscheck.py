#!/usr/bin/env python3
"""tagged_span_crosscheck.py — independent cross-checks of the tagged-span delivery.

Three checks, all read-only against committed artifacts:

  1. DENSITY CONTAINMENT. The delivered note 4.1 certifies nu_w through exact
     rational residual enclosures. We hold an independent exact rational nu_w
     from a different engine (cylinder_winding_intensity.py, wrapped by
     winding_nu_certified.py). Check that the delivered interval CONTAINS the
     independent value, and report how far the delivered centre sits from it.
     An interval that fails to contain an independently computed exact value is
     a bug in the certificate, not a rounding difference.

  2. TRUNCATION BIAS. The span-spectrum JSON measures a depth-clamped chain
     (D_MAX bins plus a tail bin), so its E[L] and CV^2 are censored. Compare it
     against the delivered all-height moments cell by cell and quantify the bias
     as a function of the reported tail fraction. This bounds how much of the
     earlier moment-limit analysis rests on censored numbers.

  3. MOMENT-LIMIT DISCRIMINATOR, UNCENSORED. Section 6.3 of the diagnostic
     conjectures Var(L)/(E L)^2 -> pi/3 - 1. Put R_w = (CV^2 - T)*w with
     T = pi/3 - 1. Model A (CV^2 = T + c/w) makes R_w a constant; model B
     (CV^2 = c/w) forces R_w to fall at slope exactly -T per unit width. The
     delivered all-height moments reach w=8 without truncation, so both models
     can be tested on clean numbers for the first time.

Run:  python scripts/tagged_span_crosscheck.py [--write]
"""
from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "geometric-consistency"

TAGGED = RESULTS / "tagged-span-resolvent.json"
SPECTRUM = RESULTS / "span-spectrum-20260913.json"
CERTIFIED = RESULTS / "winding-nu-certified-20260913.json"
OUT = RESULTS / "tagged-span-crosscheck-20260913.json"

# pi/3 - 1 to 100 significant digits, assembled as an exact rational so that
# nothing in this script depends on a floating-point constant.
PI = Fraction(
    "3.141592653589793238462643383279502884197169399375105820974944592307816406"
    "2862089986280348253421170679821480865132823066470938446095505822317253594"
)
T = PI / 3 - 1


def _centre(moments: dict) -> tuple[Fraction, Fraction, Fraction]:
    """(nu, E[L], CV^2) as exact rationals from either report shape."""
    if "centres" in moments:
        c = moments["centres"]
        return Fraction(c["nu"]), Fraction(c["mean"]), Fraction(c["cv2"])
    return Fraction(moments["nu"]), Fraction(moments["mean"]), Fraction(moments["cv2"])


def _interval(moments: dict) -> tuple[Fraction, Fraction] | None:
    iv = moments.get("density_interval")
    if not iv:
        return None
    return Fraction(iv[0]), Fraction(iv[1])


def load_tagged() -> dict[tuple[str, int, str], dict]:
    data = json.loads(TAGGED.read_text())
    out: dict[tuple[str, int, str], dict] = {}
    for system in data["systems"]:
        graph = "matching" if system["matching"] else "NN"
        for run in system["runs"]:
            key = (graph, system["width"], run["p"])
            nu, e_l, cv2 = _centre(run["moments"])
            out[key] = {
                "graph": graph,
                "width": system["width"],
                "p": run["p"],
                "tagged_states": system["tagged_states"],
                "exit_lumps": system["exit_lumps"],
                "mode": run["mode"],
                "nu": nu,
                "e_l": e_l,
                "cv2": cv2,
                "density_interval": _interval(run["moments"]),
            }
    return out


def check_density_containment(tagged: dict, certified: dict) -> dict:
    rows = []
    for (graph, w, p), rec in sorted(tagged.items()):
        mine = certified.get((graph, w, p))
        if mine is None or rec["density_interval"] is None:
            continue
        exact = Fraction(mine["nu_exact"])
        lo, hi = rec["density_interval"]
        rows.append({
            "graph": graph, "width": w, "p": p,
            "independent_nu_exact": str(exact),
            "independent_nu_float": float(exact),
            "interval_lo": str(lo), "interval_hi": str(hi),
            "interval_contains": bool(lo <= exact <= hi),
            "interval_width": float(hi - lo),
            "centre_minus_exact": float(rec["nu"] - exact),
            "centre_rel_error": float(abs(rec["nu"] - exact) / exact),
            "interval_width_rel": float((hi - lo) / exact),
        })
    return {
        "what": "delivered certified density interval vs an independent exact rational nu_w",
        "independent_engine": "scripts/cylinder_winding_intensity.py via scripts/winding_nu_certified.py",
        "points_compared": len(rows),
        "violations": sum(1 for r in rows if not r["interval_contains"]),
        "max_centre_rel_error": max((r["centre_rel_error"] for r in rows), default=None),
        "min_interval_rel_width": min((r["interval_width_rel"] for r in rows), default=None),
        "rows": rows,
    }


def check_truncation_bias(tagged: dict, spectrum: dict) -> dict:
    rows = []
    for run in spectrum["runs"]:
        graph = "matching" if run["matching"] else "NN"
        key = (graph, run["width"], run["p"])
        rec = tagged.get(key)
        if rec is None:
            continue
        dh = run["d_h_float"]
        tail = run["tail_bin_float"]
        mass = sum(dh) + tail
        e1 = sum((h + 1) * v for h, v in enumerate(dh)) / mass
        e2 = sum((h + 1) ** 2 * v for h, v in enumerate(dh)) / mass
        cv2 = e2 - e1 * e1
        cv2 /= e1 * e1
        exact_cv2 = float(rec["cv2"])
        rows.append({
            "graph": graph, "width": run["width"], "p": run["p"],
            "d_max": run["d_max"],
            "tail_fraction": tail / mass,
            "censored_e_l": e1, "all_height_e_l": float(rec["e_l"]),
            "censored_cv2": cv2, "all_height_cv2": exact_cv2,
            "cv2_rel_bias": abs(cv2 - exact_cv2) / exact_cv2,
            "R_w_censored": (cv2 - float(T)) * run["width"],
            "R_w_true": (exact_cv2 - float(T)) * run["width"],
        })
    worst = sorted(rows, key=lambda r: -r["cv2_rel_bias"])[:4]
    return {
        "what": "depth-clamped span spectrum vs delivered all-height moments, same (graph,width,p)",
        "points_compared": len(rows),
        "material_bias_threshold_tail_fraction": 1e-6,
        "material_bias_points": [r for r in rows if r["tail_fraction"] > 1e-6],
        "clean_points_max_cv2_rel_bias": max(
            (r["cv2_rel_bias"] for r in rows if r["tail_fraction"] <= 1e-6), default=None),
        "worst_four": worst,
        "all_rows": rows,
    }


def check_spectrum_vs_tagged_histogram(tagged: dict, spectrum: dict,
                                       floor_rel: float = 1e-20) -> dict:
    """Compare the two constructions height by height, not just in their moments.

    The depth-clamped chain (span_spectrum_build.cpp) tracks the ages of every
    active component and then projects onto a span histogram. The tagged
    resolvent tracks one lineage and never stores an age. If they are both
    computing the span of the same object, they must agree on d_h for every
    h <= D_MAX, and on the mass sitting above the cutoff. Comparing only the
    first two moments would not catch a compensating error; this does.

    Only heights whose value exceeds `floor_rel * nu` are scored. Below that the
    committed file stores denormal-scale float64 numbers (d_h ~ 1e-40, tail bins
    ~ 1e-45) whose relative difference is meaningless, and they are counted
    separately rather than used to characterise agreement.

    Requires importing the delivered tagged module, so it is the one check here
    that consumes delivered code rather than a committed artifact.
    """
    sys.path.insert(0, str(ROOT / "scripts"))
    import tagged_winding_span as T  # noqa: E402

    rows = []
    for run in spectrum["runs"]:
        if run["width"] > 7:
            continue  # width 8 has no committed truncated row to compare against
        graph = "matching" if run["matching"] else "NN"
        states, trans, source = T.build(run["width"], bool(run["matching"]))
        trans, source, _ = T.lump(trans, source)
        res = T.moments(trans, source, Fraction(run["p"]), bins=run["d_max"], exact=True)
        dh_exact = res["d_h"]
        dh_spec = run["d_h_float"]
        if len(dh_exact) != len(dh_spec):
            rows.append({"graph": graph, "width": run["width"], "p": run["p"],
                         "error": f"length mismatch {len(dh_exact)} vs {len(dh_spec)}"})
            continue
        nu = float(res["nu"])
        cut = floor_rel * nu
        worst_h, worst_rel = None, 0.0
        scored = below = 0
        for h, (a, b) in enumerate(zip(dh_exact, dh_spec), start=1):
            fa = float(a)
            if fa <= cut:
                below += 1
                continue
            scored += 1
            rel = abs(fa - b) / fa
            if rel > worst_rel:
                worst_h, worst_rel = h, rel
        tail_exact = float(res["tail"])
        tail_spec = run["tail_bin_float"]
        tail_scored = tail_exact > cut
        tail_rel = (abs(tail_exact - tail_spec) / tail_exact) if tail_scored else None
        rows.append({
            "graph": graph, "width": run["width"], "p": run["p"],
            "d_max": run["d_max"],
            "tagged_states": len(states), "tagged_lumps": len(trans),
            "spectrum_states": run["states"],
            "heights_total": len(dh_exact),
            "heights_scored": scored, "heights_below_floor": below,
            "floor_abs": cut,
            "max_rel_diff_d_h": (worst_rel if scored else None),
            "max_rel_diff_at_height": worst_h,
            "tail_exact": tail_exact, "tail_spectrum": tail_spec,
            "tail_scored": tail_scored, "tail_rel_diff": tail_rel,
            "d_h_tagged_head": [float(x) for x in dh_exact[:4]],
            "d_h_spectrum_head": dh_spec[:4],
        })
    scored_rows = [r for r in rows if r.get("max_rel_diff_d_h") is not None]
    worst = max((r["max_rel_diff_d_h"] for r in scored_rows), default=None)
    return {
        "what": "depth-clamped chain vs tagged resolvent, height by height",
        "constructions": {
            "spectrum": "span_spectrum_build.cpp: all component ages, then span histogram",
            "tagged": "tagged_winding_span.py: one lineage, no age, no depth cutoff",
        },
        "points_compared": len(rows),
        "points_scored": len(scored_rows),
        "scoring_floor": f"d_h > {floor_rel:g} * nu",
        "worst_rel_diff_over_all_heights": worst,
        "total_heights_below_floor": sum(r.get("heights_below_floor", 0) for r in rows),
        "state_space_ratio_note": "at w=7 the tagged lumping is 71 states against 389391; "
                                  "at w=6, 36 against 668439",
        "rows": rows,
    }


def _slope(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)


def check_model_discriminator(tagged: dict, w_min: int = 4) -> dict:
    fams: dict[tuple[str, str], dict[int, tuple[float, float]]] = {}
    for (graph, w, p), rec in tagged.items():
        fams.setdefault((graph, p), {})[w] = (float(rec["e_l"]), float(rec["cv2"]))

    tv = float(T)
    families = []
    for (graph, p), by_w in sorted(fams.items(), key=lambda kv: (kv[0][0], float(Fraction(kv[0][1])))):
        ws = sorted(w for w in by_w if w >= w_min)
        r_w = {w: (by_w[w][1] - tv) * w for w in ws}
        row = {
            "graph": graph, "p": p, "widths": ws,
            "cv2": {w: by_w[w][1] for w in ws},
            "R_w": {w: r_w[w] for w in ws},
        }
        if len(ws) >= 3:
            xs = [float(w) for w in ws]
            ys = [r_w[w] for w in ws]
            n = len(xs)
            my = sum(ys) / n
            row["slope_R_w_vs_w"] = _slope(xs, ys)
            row["model_B_required_slope"] = -tv
            row["rms_model_A_const"] = (sum((y - my) ** 2 for y in ys) / n) ** 0.5
            b0 = my + tv * (sum(xs) / n)
            row["rms_model_B_slope_minus_T"] = (
                sum((y - (b0 - tv * x)) ** 2 for x, y in zip(xs, ys)) / n) ** 0.5
            row["verdict"] = ("model A" if row["rms_model_A_const"] < row["rms_model_B_slope_minus_T"]
                              else "model B")
            row["rms_advantage"] = row["rms_model_B_slope_minus_T"] - row["rms_model_A_const"]
            row["widths_4_to_8"] = [w for w in ws if w in (4, 8)]
            if 4 in r_w and 8 in r_w:
                row["model_B_prediction_R_8"] = r_w[4] - 4 * tv
                row["measured_R_8"] = r_w[8]
                row["model_B_gap_at_w8"] = abs(row["model_B_prediction_R_8"] - r_w[8])
        families.append(row)

    testable = [f for f in families if "verdict" in f]
    return {
        "what": "section 6.3 moment-limit reading, tested on uncensored all-height moments",
        "target_T": "pi/3 - 1",
        "T_float": tv,
        "definition": "R_w = (CV^2 - T)*w ; model A -> constant, model B -> slope exactly -T",
        "families_tested": len(testable),
        "families_favouring_model_A": sum(1 for f in testable if f["verdict"] == "model A"),
        "min_model_B_gap_at_w8": min(
            (f["model_B_gap_at_w8"] for f in testable if "model_B_gap_at_w8" in f), default=None),
        "families": families,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true", help=f"write {OUT.name}")
    args = ap.parse_args()

    tagged = load_tagged()
    spectrum = json.loads(SPECTRUM.read_text())
    certified = {(r["graph"], r["width"], r["p"]): r
                 for r in json.loads(CERTIFIED.read_text())["rows"]}

    report = {
        "schema": "matching-one/tagged-span-crosscheck/1",
        "date": "2026-09-13",
        "pr": 739,
        "inputs": {
            "tagged": str(TAGGED.relative_to(ROOT)),
            "spectrum": str(SPECTRUM.relative_to(ROOT)),
            "certified": str(CERTIFIED.relative_to(ROOT)),
        },
        "density_containment": check_density_containment(tagged, certified),
        "truncation_bias": check_truncation_bias(tagged, spectrum),
        "histogram_agreement": check_spectrum_vs_tagged_histogram(tagged, spectrum),
        "model_discriminator": check_model_discriminator(tagged),
        "honesty": [
            "This script verifies arithmetic and mutual consistency. It does not "
            "verify the tagged automaton's construction, the unique-anchor pathwise "
            "theorem, or any asymptotic claim.",
            "The independent nu_w used in check 1 comes from a different exact engine "
            "but shares the same underlying cylinder model; it is an independent "
            "implementation, not an independent model.",
            "Check 3 compares two constructions of the same object, but the tagged "
            "side is delivered code imported at run time, so a shared misreading of "
            "the span definition would not be caught by it. The two state spaces are "
            "structurally unrelated, which is the reason the agreement is meaningful.",
            "Check 4 fits four widths at most. Model A beating model B on rms is a "
            "consistency statement about two one-parameter readings, not evidence "
            "for any particular correction exponent.",
        ],
    }

    d = report["density_containment"]
    print(f"[1] density containment: {d['points_compared']} points, "
          f"{d['violations']} violations, worst centre rel.err {d['max_centre_rel_error']:.2e}")
    t = report["truncation_bias"]
    print(f"[2] truncation bias: {t['points_compared']} points; "
          f"{len(t['material_bias_points'])} exceed 1e-6 tail and are materially biased; "
          f"clean points max CV^2 rel.bias {t['clean_points_max_cv2_rel_bias']:.2e}")
    h = report["histogram_agreement"]
    print(f"[3] histogram agreement: {h['points_compared']} points, "
          f"worst d_h rel.diff over all heights {h['worst_rel_diff_over_all_heights']:.2e}")
    m = report["model_discriminator"]
    print(f"[4] model discriminator: {m['families_favouring_model_A']}/{m['families_tested']} "
          f"families favour model A; smallest model-B gap at w=8 {m['min_model_B_gap_at_w8']:.4f}")

    if args.write:
        if OUT.exists():
            raise FileExistsError(f"{OUT} exists; refusing to overwrite")
        OUT.write_text(json.dumps(report, indent=1) + "\n")
        print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
