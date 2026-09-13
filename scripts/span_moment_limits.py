"""span_moment_limits.py — does Var(L)/(E L)^2 approach pi/3 - 1, or 0?

Reads the committed span-spectrum results JSON and tests the two readings of the
section-6.3 conjecture that the finite-width data can actually separate:

    model A:  Var/(E L)^2 = (pi/3 - 1) + c/w      [the conjectured constant]
    model B:  Var/(E L)^2 = c/w                   [limit 0, "ratio decays away"]

Model A fixes the limit a priori (one free parameter, c); model B also has one
free parameter. The discriminating statistic is the stability of
R*w = (Var/(E L)^2 - (pi/3 - 1)) * w: constant under model A, divergent under B.

The w=8 columns for NN p=1/4 and matching p=1/8 are the owner's one-lineage
tagged resolvent (PR #739 comment 5653178247), not this delivery's bin chain;
they are marked in the output. Run:

    python scripts/span_moment_limits.py [results.json]
"""
from __future__ import annotations
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import span_spectrum_solve as S  # noqa: E402

TARGET = math.pi / 3 - 1  # 0.04719755119659763

FAMILIES = {
    "NN p=1/4": (("NN p=1/4", "NN w=7 p=1/4"), "1/4", (8, 4.56334396378038, 0.127163500068044)),
    "NN p=1/8": (("NN p=1/8", "NN w=7 p=1/8"), "1/8", None),
    "matching p=1/8": (("matching p=1/8", "matching w=7 p=1/8"), "1/8",
                       (8, 4.711805413158092, 0.116004934830910)),
    "matching p=1/16": (("matching p=1/16", "matching w=7 p=1/16"), "1/16", None),
}


def family_series(runs, labels, p, extra_w8=None):
    rows = sorted([r for r in runs if r["label"] in labels and r["p"] == p],
                  key=lambda r: r["width"])
    ws = [r["width"] for r in rows]
    E = [S.spectrum_moments(r["d_h_float"], r["tail_bin_float"], r["d_max"])["E_L"] for r in rows]
    V = [S.spectrum_moments(r["d_h_float"], r["tail_bin_float"], r["d_max"])["var_over_E2"] for r in rows]
    if extra_w8:
        ws.append(extra_w8[0]); E.append(extra_w8[1]); V.append(extra_w8[2])
    return ws, E, V


def fit(w, V):
    """Return (c_A, rms_A, c_B, rms_B) for the two models."""
    n = len(w)
    cA = sum((v - TARGET) * wi for v, wi in zip(V, w)) / n
    rmsA = math.sqrt(sum((v - (TARGET + cA / wi)) ** 2 for v, wi in zip(V, w)) / n)
    s1 = sum(v / wi for v, wi in zip(V, w))
    s2 = sum(1.0 / (wi ** 2) for wi in w)
    cB = s1 / s2
    rmsB = math.sqrt(sum((v - cB / wi) ** 2 for v, wi in zip(V, w)) / n)
    return cA, rmsA, cB, rmsB


def slope(w, V, w_min=4):
    """Least-squares slope of R*w = (V - TARGET)*w against w, for w >= w_min.

    Discriminator: under model A (limit TARGET) this slope is 0; under model B
    (limit 0) it is exactly -TARGET per unit width, since R*w = c - TARGET*w.
    """
    pts = [(wi, (v - TARGET) * wi) for wi, v in zip(w, V) if wi >= w_min]
    n = len(pts)
    mx = sum(x for x, _ in pts) / n
    my = sum(y for _, y in pts) / n
    num = sum((x - mx) * (y - my) for x, y in pts)
    den = sum((x - mx) ** 2 for x, _ in pts)
    return num / den, n


def report(results_path=None, verbose=True):
    path = Path(results_path) if results_path else ROOT / "results" / "geometric-consistency" / "span-spectrum-20260913.json"
    with open(path) as fh:
        runs = json.load(fh)["runs"]
    out = []
    if verbose:
        print("target pi/3 - 1 = %.15f" % TARGET)
        print("model A: V = TARGET + c/w      model B: V = c/w  (limit 0)\n")
    for name, (labels, p, w8) in FAMILIES.items():
        w, E, V = family_series(runs, labels, p, w8)
        cA, rmsA, cB, rmsB = fit(w, V)
        sl, npt = slope(w, V)
        rw = [(v - TARGET) * wi for v, wi in zip(V, w)]
        late = [r for wi, r in zip(w, rw) if wi >= 4]
        spread = (max(late) - min(late)) / (sum(late) / len(late))
        src = "bin chain w<=7 + tagged w=8" if w8 else "bin chain w<=6 + w=7"
        row = {
            "family": name, "source": src,
            "w": w, "E_L": E, "var_over_E2": V,
            "R_w": rw, "R_w_spread_w_ge_4": spread,
            "model_A": {"c": cA, "rms": rmsA, "limit": TARGET},
            "model_B": {"c": cB, "rms": rmsB, "limit": 0.0},
            "slope_Rw_w_ge_4": sl, "slope_points": npt,
            "model_B_required_slope": -TARGET,
            "excludes_limit_zero": bool(sl > -0.5 * TARGET),
        }
        out.append(row)
        if not verbose:
            continue
        print("%-16s [%s]" % (name, src))
        print("   w       : " + " ".join("%7d" % wi for wi in w))
        print("   V       : " + " ".join("%7.4f" % v for v in V))
        print("   R*w     : " + " ".join("%7.4f" % r for r in rw))
        print("   [A] c=%.4f rms=%.5f | [B] c=%.4f rms=%.5f | rms ratio %.2fx | R*w spread(w>=4) %.1f%%"
              % (cA, rmsA, cB, rmsB, rmsB / rmsA, 100 * spread))
        print("   slope(R*w) for w>=%d = %+.5f  (model B requires %+.5f)  -> limit 0 %s"
              % (4, sl, -TARGET, "EXCLUDED" if row["excludes_limit_zero"] else "not excluded"))
        print()
    return out


if __name__ == "__main__":
    import json as _json
    res = report(sys.argv[1] if len(sys.argv) > 1 else None)
    if len(sys.argv) > 2:
        with open(sys.argv[2], "w") as fh:
            _json.dump({"target": TARGET, "families": res}, fh, indent=1)
