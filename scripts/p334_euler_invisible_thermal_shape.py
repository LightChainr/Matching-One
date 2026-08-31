#!/usr/bin/env python3
"""Full signed birth histograms of the frozen equal-degree safe-label tangent."""
import argparse
import csv
import hashlib
import json
import math

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import binom

from p334_safe_contact_response import ROOT, SOURCE, FORK_PATH, P_REF, blob, array_csv

CONTACT = "959a7fa26677c416b874d272f1ba66523fb38f73"
DENOMINATOR = 8 * 8000
OUT = ROOT / "results/p334-euler-invisible-thermal-shape"


def extract():
    OUT.mkdir(parents=True, exist_ok=True)
    target = OUT / "signed_birth_histograms.json"
    if target.exists():
        raise FileExistsError(target)
    result = {
        "source_commit": SOURCE, "contact_commit": CONTACT,
        "tangent_definition_commit": "9ce53a5a",
        "p_ref": P_REF, "batch_denominator": DENOMINATOR,
        "histogram_order": ["F1", "F2"],
        "estimand": "Equal orientation mixture/full original20k denominator; GROUPS[2]: own-R0, U/V own-safe and equal contact degree. Per selected quartet/orientation dg=(eU-cU)-(eV-cV), integer +dg on U0/U1 K1/K2 and -dg on V0/V1, batch divisor64000.",
        "inference": "Original20 batches, not independent thresholds/tails/orientations.",
        "source_sha256": {}, "sizes": {},
    }
    for n in (325, 425):
        all_hist, counts = [], []
        for batch in range(20):
            raw_path = f"{FORK_PATH}/N{n}/N{n}.batch{batch:02}.csv.gz"
            mark_path = f"results/p334-next-label-contact-coordinates/N{n}/N{n}.batch{batch:02}.csv.gz"
            raw_blob, mark_blob = blob(SOURCE, raw_path), blob(CONTACT, mark_path)
            h, raw = array_csv(raw_blob)
            ch, marks = array_csv(mark_blob)
            hi, ci = {x: i for i, x in enumerate(h)}, {x: i for i, x in enumerate(ch)}
            raw = raw[np.lexsort(tuple(raw[:, hi[k]] for k in ("replica", "group", "quartet", "counter")))].reshape(1000, 8, 2, 2, -1)
            marks = marks[np.lexsort(tuple(marks[:, ci[k]] for k in ("group", "quartet", "counter")))].reshape(1000, 8, 2, -1)
            for key in ("counter", "quartet", "group", "next_label"):
                if not np.array_equal(raw[..., 0, hi[key]], marks[..., ci[key]]):
                    raise ValueError(f"Source alignment N{n} batch{batch}: {key}")
            hist = np.zeros((2, n + 1), dtype=np.int64)
            batch_counts = []
            for orientation in ("first", "second"):
                old = raw[:, 0, 0, 0, hi[f"{orientation}_rank"]]
                e, c = (marks[..., ci[f"{orientation}_{key}"]] for key in ("e", "c"))
                nr = raw[..., 0, hi[f"{orientation}_next_rank"]]
                chosen = ((old == 0)[:, None] & (nr[:, :, 0] == 0)
                          & (nr[:, :, 1] == 0) & (e[:, :, 0] == e[:, :, 1]))
                dg = ((e-c)[:, :, 0]-(e-c)[:, :, 1])[chosen]
                batch_counts.append({"selected_quartets": int(chosen.sum()),
                                     "nonzero_dg_quartets": int(np.count_nonzero(dg))})
                for j, key in enumerate(("k1", "k2")):
                    ks = raw[..., hi[f"{orientation}_{key}"]][chosen]
                    for group, sign in ((0, 1), (1, -1)):
                        for replica in (0, 1):
                            np.add.at(hist[j], ks[:, group, replica], sign * dg)
            all_hist.append(hist.tolist())
            counts.append(batch_counts)
            result["source_sha256"][raw_path] = hashlib.sha256(raw_blob).hexdigest()
            result["source_sha256"][mark_path] = hashlib.sha256(mark_blob).hexdigest()
        result["sizes"][str(n)] = {
            "batch_ids": list(range(20)), "batch_integer_histograms": all_hist,
            "selected_counts_first_second": counts,
        }
        target.write_text(json.dumps(result, indent=2) + "\n")
        print(f"N{n}: full20batch integer F1/F2 histograms saved", flush=True)


def reduced_roots(integer_coefficients):
    """Remove exact endpoint factors before ordinary double-precision root readout.

    This avoids declaring underflow/roundoff zeros to be thermal roots. It is
    not a high-precision all-roots certificate; no new curve model is fitted.
    """
    b = np.asarray(integer_coefficients)
    ix = np.flatnonzero(b)
    if not len(ix):
        return {"identically_zero": True, "roots": []}
    n, lo, hi = len(b)-1, int(ix[0]), int(ix[-1])
    degree = hi-lo
    q = np.array([float(int(b[lo+j])*math.comb(n, lo+j)/math.comb(degree, j))
                  for j in range(degree+1)])
    q /= np.max(np.abs(q))
    kk = np.arange(degree+1)
    evaluate = lambda p: float(binom.pmf(kk, degree, p) @ q)
    grid = np.linspace(0, 1, 2001)
    values = binom.pmf(kk[None, :], degree, grid[:, None]) @ q
    intervals = [(float(grid[i]), float(grid[i+1])) for i in range(len(grid)-1)
                 if np.sign(values[i]) * np.sign(values[i+1]) < 0]
    roots = [float(brentq(evaluate, a, b, xtol=2e-14)) for a, b in intervals]
    return {"identically_zero": False, "endpoint_factors": [lo, n-hi],
            "reduced_degree": degree, "roots": roots,
            "method": "Exact zero-endpoint coefficients factored; fixed2001-point reduced-Bernstein sign brackets + Brent. Ordinary numerical crossings, not certified absence of tangent/even-multiplicity roots."}


def score():
    source_path = OUT / "signed_birth_histograms.json"
    raw = json.loads(source_path.read_text())
    result = {"histogram_file": str(source_path.relative_to(ROOT)),
              "histogram_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
              "labels": ["F1", "F2", "A", "E"], "p_ref": P_REF,
              "definition": "A=F1+F2-1; E=1-F1+F2. Tangent constants vanish. Curves use exact integer cumulative histograms as Bernstein coefficients before float normalization.",
              "covariance": "Every readout uses the same original20batch vectors. The saved integerhistograms reconstruct the complete joint curve covariance; no p-point/tail independence and no simultaneous confidence-band claim.",
              "sizes": {}}
    pgrid = np.unique(np.r_[np.linspace(0, 1, 1001), P_REF])
    for ns, src in raw["sizes"].items():
        n = int(ns)
        hist = np.asarray(src["batch_integer_histograms"], dtype=np.int64)
        c12 = np.cumsum(hist, axis=2)  # Exact zero tails before any division.
        coef = np.stack((c12[:, 0], c12[:, 1], c12[:, 0]+c12[:, 1],
                         c12[:, 1]-c12[:, 0]), axis=1)
        kk = np.arange(n+1)

        def values_at(p):
            return (coef @ binom.pmf(kk, n, p))/DENOMINATOR

        def readout_at(p):
            b = values_at(p)
            return {"p": float(p), "mean": b.mean(0).tolist(),
                    "se": (b.std(0, ddof=1)/np.sqrt(20)).tolist(),
                    "joint20_batch_values": b.tolist()}

        ec = coef[:, 3].sum(0)
        root_info = reduced_roots(ec)
        roots = root_info["roots"]
        main = next(r for r in roots if r > P_REF)
        ecurve = lambda p: float(values_at(p)[:, 3].mean())
        derivative = n*np.diff(coef[:, 3], axis=1)/DENOMINATOR
        slope = float((derivative @ binom.pmf(np.arange(n), n-1, main)).mean())
        root_readout = readout_at(main)
        local_se = root_readout["se"][3]/abs(slope)
        # A fixed local interval reports missing/ambiguous leave-one-batch roots
        # rather than moving to a compatible branch.
        loo_sets = [reduced_roots(ec-row)["roots"] for row in coef[:, 3]]
        loo_main = [[r for r in rr if .55 < r < .70] for rr in loo_sets]
        valid = [rr[0] for rr in loo_main if len(rr) == 1]
        loo_se = float(np.sqrt(19/20*np.sum((np.array(valid)-np.mean(valid))**2))) if len(valid) == 20 else None
        main_data = {"p": main, "slope": slope, "local_delta_se": local_se,
                     "LOO_root_sets": loo_sets, "LOO_main_window": [.55, .70],
                     "LOO_main_roots": [rr[0] if len(rr) == 1 else None for rr in loo_main],
                     "LOO_main_valid_count": len(valid), "LOO_se": loo_se}
        bounds = [0., *roots, 1.]
        lobes = []
        for left, right in zip(bounds[:-1], bounds[1:]):
            sign = 1 if ecurve((left+right)/2) > 0 else -1
            opt = minimize_scalar(lambda p: -sign*ecurve(p), bounds=(left, right),
                                  method="bounded", options={"xatol": 2e-12})
            point = readout_at(opt.x)
            # Integral Bernstein basis = beta CDF/(N+1), evaluated jointly.
            from scipy.special import betainc
            integrated_basis = (betainc(kk+1, n-kk+1, right)-betainc(kk+1, n-kk+1, left))/(n+1)
            batch_area = (coef[:, 3] @ integrated_basis)/DENOMINATOR
            lobes.append({"interval": [left, right], "sign": sign,
                          "extremum": point, "signed_area": float(batch_area.mean()),
                          "area_se_at_point_root_limits": float(batch_area.std(ddof=1)/np.sqrt(20))})
        basis = binom.pmf(kk[:, None], n, pgrid[None, :])
        curves = (coef @ basis)/DENOMINATOR
        means, ses = curves.mean(0), curves.std(0, ddof=1)/np.sqrt(20)
        clock12 = (hist @ kk)/DENOMINATOR
        integrals = coef.sum(2)/(DENOMINATOR*(n+1))
        key_readouts = {"p_ref": readout_at(P_REF), "at_main_root": root_readout,
                        "integrals": {"mean": integrals.mean(0).tolist(),
                                      "se": (integrals.std(0, ddof=1)/np.sqrt(20)).tolist(),
                                      "joint20_batch_values": integrals.tolist()},
                        "K1_K2": {"mean": clock12.mean(0).tolist(),
                                  "se": (clock12.std(0, ddof=1)/np.sqrt(20)).tolist()}}
        result["sizes"][ns] = {"batch_ids": src["batch_ids"], "root_readout": root_info,
                               "main_root": main_data, "lobes": lobes,
                               "key_readouts": key_readouts,
                               "curve_grid_p": pgrid.tolist(), "curve_mean": means.tolist(),
                               "curve_pointwise_se": ses.tolist()}
        print(f"N{n}: roots={roots}; main={main:.9f}+/-{local_se:.8f} local, LOO={loo_se}", flush=True)
        for lobe in lobes:
            v = lobe["extremum"]
            print(f"  sign{lobe['sign']:+}: p={v['p']:.8f} E={v['mean'][3]:.9g}+/-{v['se'][3]:.6g}; area={lobe['signed_area']:.9g}", flush=True)
    (OUT/"score.json").write_text(json.dumps(result, indent=2)+"\n")
    with (OUT/"thermal_curve.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["N", "p", "prefixes", "original_batches",
                         *[f"mean_{s}" for s in result["labels"]],
                         *[f"pointwise_se_{s}" for s in result["labels"]]])
        for n, row in result["sizes"].items():
            for i, p in enumerate(row["curve_grid_p"]):
                writer.writerow([n, p, 20000, 20,
                                 *[a[i] for a in row["curve_mean"]],
                                 *[a[i] for a in row["curve_pointwise_se"]]])


def plot():
    """Repository-export figure: shape and common pointwise batch uncertainty."""
    import matplotlib as mpl
    mpl.use("Agg")
    mpl.rcParams["svg.fonttype"] = "none"
    import matplotlib.pyplot as plt
    r = json.loads((OUT/"score.json").read_text())
    contract = {
        "question": "Where does the fixed Euler-invisible tangent change the two-birth response?",
        "takeaway": "Early positive and stronger later negative E lobes arise from a switch in which delayed birth dominates.",
        "family": "line, two-size small multiples with pointwise uncertainty",
        "data": "1002 p evaluations of exact integer-histogram Bernstein polynomials;20originalbatch vectors perN, equal orientations/full20kprefix denominator",
        "surface": "Repository static PNG/SVG export, generated from score.json",
        "palette": {"policy": "hard two-root cap", "blue": "#2563A6", "orange": "#C87823", "neutral": "#444444"},
        "noncolor": "F1 dashed/F2 solid; zero and p_ref reference; separate size panels",
        "footprint": "10x6.8inch; focused p .52-.82, all-p numerical data retained",
        "uncertainty": "One SE pointwise at fixedp, same20batch covariance, not simultaneous bands or peak-search intervals",
    }
    (OUT/"chart_contract.json").write_text(json.dumps(contract, indent=2)+"\n")
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.8), sharex=True, sharey="row")
    blue, orange, neutral = "#2563A6", "#C87823", "#444444"
    for col, (ns, x) in enumerate(r["sizes"].items()):
        p = np.array(x["curve_grid_p"])
        means, ses = np.array(x["curve_mean"])*1e4, np.array(x["curve_pointwise_se"])*1e4
        axes[0, col].plot(p, means[3], color=blue, lw=2, label=r"$\partial_t E(p)$")
        axes[0, col].fill_between(p, means[3]-ses[3], means[3]+ses[3], color=blue, alpha=.16)
        root = x["main_root"]["p"]
        axes[0, col].axvline(root, color=neutral, ls=":", lw=.9)
        axes[0, col].set_title(f"N{ns}  |  main crossing {root:.4f}", loc="left", fontsize=11)
        for i, color, style in ((0, orange, "--"), (1, blue, "-")):
            axes[1, col].plot(p, means[i], color=color, ls=style, lw=1.8,
                              label=rf"$\partial_t F_{i+1}(p)$")
            axes[1, col].fill_between(p, means[i]-ses[i], means[i]+ses[i], color=color, alpha=.10)
        for row in (0, 1):
            ax = axes[row, col]
            ax.axhline(0, color=neutral, lw=.8)
            ax.axvline(P_REF, color=neutral, ls="--", lw=.8, alpha=.65)
            ax.set_xlim(.52, .82)
            ax.spines[["top", "right"]].set_visible(False)
            ax.grid(axis="y", color="#DDDDDD", lw=.5)
            ax.legend(frameon=False, fontsize=9, loc="lower right")
        axes[1, col].set_xlabel("Occupation probability p")
    axes[0, 0].set_ylabel(r"Even-topology tangent ($10^{-4}$)")
    axes[1, 0].set_ylabel(r"Birth CDF tangents ($10^{-4}$)")
    fig.suptitle("Euler-invisible next-label thermal response", x=.08, ha="left", fontsize=15)
    fig.text(.08, .935, "Equal orientation mixture; full20k prefixes perN; original20batch pointwise ±1 SE", fontsize=9, color=neutral)
    fig.text(.975, .975, "✿", ha="right", va="top", color=blue, fontsize=13)
    fig.text(.08, .02, "Dashed vertical: p_ref. Dotted: main E crossing. All-p curves and signed integer histograms are archived.", fontsize=8.5, color=neutral)
    fig.tight_layout(rect=(.015, .05, .995, .915))
    for suffix in ("png", "svg"):
        fig.savefig(OUT/f"thermal_shape.{suffix}", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["extract", "score", "plot"])
    args = parser.parse_args()
    {"extract": extract, "score": score, "plot": plot}[args.mode]()


if __name__ == "__main__":
    main()
