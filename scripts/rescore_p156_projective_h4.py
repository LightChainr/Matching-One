#!/usr/bin/env python3
"""Zero-new-sampling rescore of the archived #156 N30/N56 primitive pilot.

The archive retains batch counts for the three named primitive lines l0,l1,l2
but aggregates every other rank-one line into rank1_other.  Therefore we can
recover exactly, with the original batch covariance:

  * the frozen real C3 character C_C3;
  * the physically embedded three-line spin-4 projection A4_three(tau).

We CANNOT reconstruct the full all-primitive A4_full from this archive because
the directions inside rank1_other were not retained.  The script reports this
explicitly and never imputes those missing phases.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from statistics import mean

import mpmath as mp

from pinson_arguin_primitive import primitive_probability_direct


DESIGNS = {
    "pell_Dminus2_N30": mp.mpc(mp.mpf("0.5"), mp.mpf(5) / 6),
    "pell_Dplus1_N56": mp.mpc(mp.mpf("0.5"), mp.mpf(7) / 8),
}


def _complex_stats(values: list[complex]) -> dict[str, float]:
    n = len(values)
    if n < 2:
        raise ValueError("need at least two independent batches")
    mr = sum(z.real for z in values) / n
    mi = sum(z.imag for z in values) / n
    vr = sum((z.real - mr) ** 2 for z in values) / (n - 1)
    vi = sum((z.imag - mi) ** 2 for z in values) / (n - 1)
    cov = sum((z.real - mr) * (z.imag - mi) for z in values) / (n - 1)
    return {
        "mean_real": mr,
        "mean_imag": mi,
        "se_real": math.sqrt(vr / n),
        "se_imag": math.sqrt(vi / n),
        "cov_mean_real_imag": cov / n,
        "z_real": mr / math.sqrt(vr / n) if vr else float("inf"),
        "z_imag": mi / math.sqrt(vi / n) if vi else float("inf"),
    }


def _real_stats(values: list[float]) -> dict[str, float]:
    n = len(values)
    m = mean(values)
    variance = sum((x - m) ** 2 for x in values) / (n - 1)
    se = math.sqrt(variance / n)
    return {"mean": m, "se": se, "z": m / se if se else float("inf")}


def continuum_three_line(tau: mp.mpc) -> dict[str, object]:
    # engine l0=(1,0), l1=(0,1), l2=(1,-1)
    # maps to paper sectors {1,0}, {0,-1}, {1,1}.
    p0 = primitive_probability_direct(1, 0, tau)
    p1 = primitive_probability_direct(0, -1, tau)
    p2 = primitive_probability_direct(1, 1, tau)
    z0 = mp.mpc(1)
    z1 = (tau / abs(tau)) ** 4
    z2 = ((1 - tau) / abs(1 - tau)) ** 4
    a4 = p0 * z0 + p1 * z1 + p2 * z2
    c3 = p0 - (p1 + p2) / 2
    q3 = mp.sqrt(3) * (p2 - p1) / 2
    return {
        "probabilities": [str(p0), str(p1), str(p2)],
        "z4_weights": [[float(mp.re(z0)), float(mp.im(z0))],
                       [float(mp.re(z1)), float(mp.im(z1))],
                       [float(mp.re(z2)), float(mp.im(z2))]],
        "A4_three_real": float(mp.re(a4)),
        "A4_three_imag": float(mp.im(a4)),
        "C_C3": float(c3),
        "Q_C3": float(q3),
    }


def rescore(csv_path: Path, dps: int) -> dict[str, object]:
    mp.mp.dps = dps
    grouped: dict[str, list[dict[str, int]]] = {name: [] for name in DESIGNS}
    with csv_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            design = row["design"]
            if design not in grouped:
                continue
            grouped[design].append({
                key: int(row[key])
                for key in ("samples", "l0", "l1", "l2", "rank1_other")
            })

    out: dict[str, object] = {
        "source": str(csv_path),
        "missing_information": (
            "rank1_other is aggregated, so A4_full over all primitive slopes "
            "is not recoverable from this archive"
        ),
        "designs": {},
    }

    for design, tau in DESIGNS.items():
        rows = grouped[design]
        if not rows:
            raise ValueError(f"no batches found for {design}")
        continuum = continuum_three_line(tau)
        z0 = 1 + 0j
        z1 = complex((tau / abs(tau)) ** 4)
        z2 = complex(((1 - tau) / abs(1 - tau)) ** 4)

        a4_residual_batches: list[complex] = []
        c3_residual_batches: list[float] = []
        other_fractions: list[float] = []
        for row in rows:
            n = row["samples"]
            empirical_a4 = (
                row["l0"] * z0 + row["l1"] * z1 + row["l2"] * z2
            ) / n
            continuum_a4 = complex(
                continuum["A4_three_real"], continuum["A4_three_imag"]
            )
            a4_residual_batches.append(empirical_a4 - continuum_a4)

            empirical_c3 = (row["l0"] - (row["l1"] + row["l2"]) / 2) / n
            c3_residual_batches.append(empirical_c3 - continuum["C_C3"])
            other_fractions.append(row["rank1_other"] / n)

        out["designs"][design] = {
            "tau": [float(mp.re(tau)), float(mp.im(tau))],
            "batches": len(rows),
            "samples_per_batch": rows[0]["samples"],
            "continuum_three_line": continuum,
            "A4_three_residual": _complex_stats(a4_residual_batches),
            "C_C3_residual": _real_stats(c3_residual_batches),
            "rank1_other_fraction": _real_stats(other_fractions),
        }

    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("results/local-20260829/P156-square-bond-primitive-pilot/result.batches.csv"),
    )
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = rescore(args.input, args.dps)
    text = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)
