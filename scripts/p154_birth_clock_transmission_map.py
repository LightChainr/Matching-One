#!/usr/bin/env python3
"""Evaluate an analytic birth-clock map from already saved unmarked p-jets.

No source response is fitted, no permutations are replayed, and no new
prospective data are read. Output gains are nominal planning values.
"""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

BASE_REV = "7da1eeb0"
BASE_PATH = "results/norm4-source-endpoint-1m/latest.json"
MODES = {
    "00_common": ([1., 1.], [1., 1.]),
    "10_orientation_common_birth": ([1., -1.], [1., -1.]),
    "01_relative_birth": ([-1., -1.], [1., 1.]),
    "11_orientation_relative_birth": ([-1., 1.], [1., -1.]),
}


def evaluate(n, baseline):
    rows = [baseline["direction"][g] for g in ("first", "second")]
    q, e, qq, ee = [np.array([r[k] for r in rows])
                     for k in ("q_p", "E_p", "q_pp", "E_pp")]
    d, t = q.mean(), qq.mean()
    a = n ** (13 / 8) / 2
    delta = baseline["delta_cos4"]
    p4 = lambda v: (v[0] - v[1]) / delta
    f = ((q - e) / 2, (q + e) / 2)
    fp = ((qq - ee) / 2, (qq + ee) / 2)
    modes = {}
    for name, alpha in MODES.items():
        alpha = tuple(np.array(v) for v in alpha)
        gains = {}
        for tag, value, slope in (("value", 1., 0.), ("p_derivative", 0., 1.)):
            j = tuple(value * al * ff / n for al, ff in zip(alpha, f))
            jp = tuple(al * (value * ffp + slope * ff) / n
                       for al, ff, ffp in zip(alpha, f, fp))
            jq, jqp = j[0] + j[1], jp[0] + jp[1]
            rootdot = -jq.mean() / d
            ddot = jqp.mean() + rootdot * t
            channels = [float(sign * a / d *
                              (p4(jjp) + rootdot * p4(ffp) - p4(ff) * ddot / d))
                        for sign, ff, ffp, jjp in zip((-1, 1), f, fp, jp)]
            gains[tag] = dict(zip(("entry", "completion"), channels))
            gains[tag]["net"] = sum(channels)
        cv = gains["value"]
        gains["flat_mode_completion_over_entry"] = (
            cv["completion"] / cv["entry"] if name != "00_common" else None)
        gains["flat_mode_amplitude_for_net_0p5"] = (
            .5 / cv["net"] if name != "00_common" else None)
        modes[name] = gains
    return {"N": n, "p0": baseline["p0"], "D": float(d),
            "delta_cos4": delta, "baseline_p_jets": {
                g: {k: baseline["direction"][g][k]
                    for k in ("q_p", "E_p", "q_pp", "E_pp")}
                for g in ("first", "second")}, "gains": modes}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    commit = subprocess.check_output(["git", "rev-parse", BASE_REV], text=True).strip()
    raw = subprocess.check_output(["git", "show", commit + ":" + BASE_PATH])
    source = json.loads(raw)
    output = {
        "status": "ANALYTIC_MAP_WITH_NOMINAL_UNMARKED_BASELINE_GAINS",
        "baseline_commit": commit, "baseline_path": BASE_PATH,
        "baseline_sha256": hashlib.sha256(raw).hexdigest(),
        "source_coefficients_fitted": False, "new_samples": 0,
        "baseline_uncertainty_propagated": False,
        "scope": "Exact linear map; numeric gains are old-baseline point calibrations, not prospective intervals or source detections.",
        "alpha_definition": "alpha_jg=N*J_jg/Fprime_jg=N*p*E_birth_posterior[z_(k-1)/k]",
        "mode_definition": "alpha_jg=a00+eta_g*a10+tau_j*a01+eta_g*tau_j*a11; eta=(+1,-1),tau=(-1,+1)",
        "by_N": {str(n): evaluate(n, source["by_N"][str(n)]["source"])
                 for n in (85, 260, 340)},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(output, stream, indent=2, allow_nan=False)
        stream.write("\n")
    for n, row in output["by_N"].items():
        print(n, {k: (v["value"]["net"], v["p_derivative"]["net"],
                       v["flat_mode_completion_over_entry"])
                  for k, v in row["gains"].items()})


if __name__ == "__main__":
    main()
