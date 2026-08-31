#!/usr/bin/env python3
"""One original-batch join of Hessians, saved corners and saved new64 first means."""
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OLD = "172fbeb1ed28019b04f132859da3ea247942bb1c"
OLD_DIR = "results/p334-prefix-response-projection-joint"
HESSIAN = "c48fa360a37a9887ef32ff6d3ce947c4e4601b53"
HESSIAN_PATH = "results/p334-mixed-source-curvature/score.json"
RECTANGLE = "e747323340ea1e672c87d2ec46135ff8a75d1a08"
RECTANGLE_PATH = "results/p334-mixed-source-rectangle/score.json"
NEW = "8ad30617b0a3076a5c01a208eb213096d8879b32"
FIRST_PATH = "experiments/p334-mechanism-response-20260831/results-extension/prefix_statistics_N{n}.npz"
OUT = ROOT / "results/p334-mixed-source-curvature-joint"
ORIS = ("first", "second")
OBS = ("A_ref", "E_ref", "C", "W")
FIRST_OBS = {"A_ref": ("p_ref.A", 1.), "E_ref": ("p_ref.E", 1.),
             "C": ("p_integral.A", -.5), "W": ("p_integral.E", -1.)}
DELTA = {"325": -.7634556213017751, "425": -.8928996539792388}


def factor(loo):
    return np.sqrt(19 / 20) * (loo - loo.mean(axis=0))


def main():
    hashes = {}

    def read(commit, path):
        b = subprocess.check_output(["git", "show", f"{commit}:{path}"], cwd=ROOT)
        hashes[f"{commit}:{path}"] = sha256(b).hexdigest()
        return b

    old = json.loads(read(OLD, OLD_DIR + "/score.json"))
    hs = json.loads(read(HESSIAN, HESSIAN_PATH))
    rs = json.loads(read(RECTANGLE, RECTANGLE_PATH))
    if hs["new64_commit"] != NEW or rs["source_commit"] != NEW:
        raise ValueError("Supplied tensors/corners do not refer to the frozen new64 archive")
    OUT.mkdir(parents=True, exist_ok=True)
    sizes = {}
    lines = ["# Mixed-source curvature: one original-20-batch covariance join", "",
             "All ± values are paired delete-one-original-batch SE. No new sampling,",
             "raw fork replay, weight calculation or model fitting is performed.", ""]
    for n in ("325", "425"):
        o, h, r = old["sizes"][n], hs["sizes"][n], rs["sizes"][n]
        if not (o["batch_ids"] == h["batch_ids"] == r["batch_ids"] == list(range(20))):
            raise ValueError("Original batch deletion rows are not aligned")
        if o["prefix_counts"] != [1000] * 20 or r["population_per_batch"] != 1000:
            raise ValueError("The full original-population denominator changed")
        saved = json.loads(gzip.decompress(read(OLD, OLD_DIR + "/" + o["complete_covariance_factor_file"])))
        hp, hl = np.asarray(h["estimate"]), np.asarray(h["LOO"])
        rp, rl = np.asarray(r["estimate"]), np.asarray(r["LOO"])
        hi, ri = ({k: i for i, k in enumerate(z["labels"])} for z in (h, r))

        # Read only the saved first-response means; no determinant or raw-tail readout.
        with np.load(io.BytesIO(read(NEW, FIRST_PATH.format(n=n))), allow_pickle=False) as z:
            batch = z["batch"]
            source_labels = list(z["labels"])
            matrix = z["new64"]
            first_columns, first_labels = [], []
            for receiver in ORIS:
                for source in ORIS:
                    for obs in OBS:
                        endpoint, scale = FIRST_OBS[obs]
                        k = source_labels.index(f"{endpoint}.mean_J[{receiver},{source}]")
                        first_labels.append(f"new64.00.{receiver}.{obs}.H_{source}")
                        first_columns.append([scale * matrix[batch == b, k].sum() / 1000 for b in range(20)])
            counts = [int(np.count_nonzero(batch == b)) for b in range(20)]
        first_batch = np.column_stack(first_columns)
        fp = first_batch.mean(axis=0)
        fl = (20 * fp - first_batch) / 19
        fi = {k: i for i, k in enumerate(first_labels)}

        names, values, loos = [], [], []

        def add(name, point, loo):
            names.append(name); values.append(float(point)); loos.append(loo)

        def hc(stream, receiver, obs, comp, group="00"):
            k = hi[f"{stream}.{group}.{receiver}.{obs}.{comp}"]
            return hp[k], hl[:, k]

        def mixed(stream, receiver, obs):
            if receiver in ORIS:
                return hc(stream, receiver, obs, "fs")
            f, s = (hc(stream, ori, obs, "fs") for ori in ORIS)
            weights = (.5, .5) if receiver == "S" else (1 / DELTA[n], -1 / DELTA[n])
            return weights[0] * f[0] + weights[1] * s[0], weights[0] * f[1] + weights[1] * s[1]

        # Fixed physically named pure-versus-cross contrasts, only original00/new64.
        for receiver in ORIS:
            own, other = ("ff", "ss") if receiver == "first" else ("ss", "ff")
            for obs in OBS:
                op, ol = hc("new64", receiver, obs, own)
                for label, component in (("other_pure", other), ("mixed", "fs")):
                    cp, cl = hc("new64", receiver, obs, component)
                    add(f"new64.00.{receiver}.{obs}.own_minus_{label}", op - cp, ol - cl)

        # Mixed scores vanish off00. The finite rectangle and both Hessians keep
        # the same full1000-per-batch normalization, rather than conditioning by mass.
        for receiver in (*ORIS, "S", "D"):
            for obs in OBS:
                oldp, oldl = mixed("old8", receiver, obs)
                newp, newl = mixed("new64", receiver, obs)
                k = ri[f"all.{receiver}.{obs}.mixed_rectangle"]
                add(f"old8.{receiver}.{obs}.finite_rectangle_minus_Hfs0", rp[k] - oldp, rl[:, k] - oldl)
                add(f"{receiver}.{obs}.new64_minus_old8_Hfs0", newp - oldp, newl - oldl)

        ratio_support = {}
        for receiver in ORIS:
            own = "ff" if receiver == "first" else "ss"
            for obs in ("A_ref", "C"):
                num, numl = hc("new64", receiver, obs, own)
                k = fi[f"new64.00.{receiver}.{obs}.H_{receiver}"]
                den, denl = fp[k], fl[:, k]
                if den == 0 or np.any(denl * den <= 0):
                    raise ValueError("A preselected own-response denominator is unresolved/sign-changing")
                ratio_support[f"{receiver}.{obs}"] = {"H1": float(den), "H1_se": float(np.linalg.norm(factor(fl)[:, k])),
                    "H1_LOO_range": [float(denl.min()), float(denl.max())], "H2": float(num),
                    "H2_se": float(np.linalg.norm(factor(numl[:, None]))) }
                add(f"new64.00.{receiver}.{obs}.own_H2_over_H1", num / den, numl / denl)
                add(f"new64.00.{receiver}.{obs}.t_half_quadratic_over_linear", num / (4 * den), numl / (4 * denl))
        derived_loo = np.column_stack(loos)
        derived_factor = factor(derived_loo)

        # Crucial sign: h.factor is raw-batch centered; use its supplied negative
        # LOO_factor. Rectangle.factor and every predecessor are already LOO-centered.
        hfactor = np.asarray(h["LOO_factor"])
        rfactor = np.asarray(r["factor"])
        blocks = [np.asarray(saved["factor"]), hfactor, rfactor, factor(fl), derived_factor]
        joint = np.column_stack(blocks)
        full_labels = (saved["labels"] + ["mixed_hessian." + k for k in h["labels"]]
                       + ["mixed_rectangle." + k for k in r["labels"]]
                       + ["mixed_first_response." + k for k in first_labels]
                       + ["mixed_joint." + k for k in names])
        offsets = np.cumsum([0] + [x.shape[1] for x in blocks]).tolist()
        packed = {"batch_ids": list(range(20)), "labels": full_labels, "factor": joint.tolist(),
                  "convention": "sqrt(19/20)*(LOO-mean_LOO) throughout; covariance=factor.T@factor; rank<=19; no inverse",
                  "blocks": {k: offsets[i:i + 2] for i, k in enumerate(("previous", "hessian", "rectangle", "new64_first_means", "derived"))},
                  "hessian_labels": h["labels"], "hessian_estimate": h["estimate"], "hessian_LOO": h["LOO"],
                  "rectangle_labels": r["labels"], "rectangle_estimate": r["estimate"], "rectangle_LOO": r["LOO"],
                  "first_labels": first_labels, "first_20_batch_means": first_batch.tolist(), "first_estimate": fp.tolist(), "first_LOO": fl.tolist(),
                  "derived_labels": names, "derived_estimate": values, "derived_LOO": derived_loo.tolist()}
        filename = f"N{n}.complete_common_factor.json.gz"
        b = gzip.compress((json.dumps(packed, separators=(",", ":"), allow_nan=False) + "\n").encode(), mtime=0)
        (OUT / filename).write_bytes(b)
        sizes[n] = {"batch_ids": list(range(20)), "original_prefix_counts": [1000] * 20,
                    "new64_00_prefix_counts": counts, "labels": names, "estimate": values,
                    "LOO": derived_loo.tolist(), "factor": derived_factor.tolist(),
                    "se": np.linalg.norm(derived_factor, axis=0).tolist(),
                    "covariance": (derived_factor.T @ derived_factor).tolist(), "ratio_denominators": ratio_support,
                    "hessian_labels": h["labels"], "hessian_estimate": h["estimate"], "hessian_se": h["se"],
                    "first_labels": first_labels, "first_estimate": fp.tolist(), "first_se": np.linalg.norm(factor(fl), axis=0).tolist(),
                    "complete_covariance_factor_file": filename, "complete_covariance_factor_sha256": sha256(b).hexdigest(),
                    "complete_coordinate_count": len(full_labels)}
        lines += [f"## N{n}", "", "### New64 original00 tensor", "",
                  "| Receiver / observable | Own pure | Other pure | Mixed | Own − other | Own − mixed |",
                  "|---|---:|---:|---:|---:|---:|"]
        ni = {k: i for i, k in enumerate(names)}

        def ht(receiver, obs, comp):
            k = hi[f"new64.00.{receiver}.{obs}.{comp}"]
            return f"{hp[k]:.7g} ± {np.linalg.norm(hfactor[:, k]):.4g}"

        def dt(key):
            k = ni[key]
            return f"{values[k]:.7g} ± {np.linalg.norm(derived_factor[:, k]):.4g}"

        for receiver in ORIS:
            own, other = ("ff", "ss") if receiver == "first" else ("ss", "ff")
            for obs in OBS:
                cols = [ht(receiver, obs, c) for c in (own, other, "fs")]
                cols += [dt(f"new64.00.{receiver}.{obs}.own_minus_{c}") for c in ("other_pure", "mixed")]
                lines.append(f"| {receiver} / {obs} | " + " | ".join(cols) + " |")
        lines += ["", "### Paired finite-vs-zero-source and conditional-stream differences", "",
                  "| Observer / observable | old8 rectangle − old8 Hfs(0) | new64 Hfs(0) − old8 Hfs(0) |",
                  "|---|---:|---:|"]
        for receiver in (*ORIS, "S", "D"):
            for obs in OBS:
                lines.append(f"| {receiver} / {obs} | " + dt(f"old8.{receiver}.{obs}.finite_rectangle_minus_Hfs0")
                             + " | " + dt(f"{receiver}.{obs}.new64_minus_old8_Hfs0") + " |")
        lines += ["", "### Own-source local gain curvature", "",
                  "| Receiver / observable | H1 | H2/H1 | t=1/2 quadratic / linear |", "|---|---:|---:|---:|"]
        for receiver in ORIS:
            for obs in ("A_ref", "C"):
                q = ratio_support[f"{receiver}.{obs}"]
                lines.append(f"| {receiver} / {obs} | {q['H1']:.7g} ± {q['H1_se']:.4g} | "
                             + dt(f"new64.00.{receiver}.{obs}.own_H2_over_H1") + " | "
                             + dt(f"new64.00.{receiver}.{obs}.t_half_quadratic_over_linear") + " |")
        lines.append("")
    boundary = (
        "Fixed commuting physical-source coordinates and normalization. Pure/mixed contrasts and ratios use new64 original00, "
        "zero padded to the original population; old8 finite rectangle minus old8 Hfs is paired within its original stream. "
        "New64-minus-old8 keeps shared prefixes/batches and is not independent population replication. The side-one rectangle "
        "is an integral of Hfs over [-1/2,1/2]^2, not identically Hfs(0). H2/H1 is a local fractional gain derivative; "
        "one-quarter of it is the ratio of specified second/first Taylor contributions at t=1/2, not a finite-response prediction. "
        "No ratios are formed for cross responses, E, or W. All old/new matrices and first means are appended jointly to172fbeb1 "
        "with consistently signed original20 LOO factors. No fit, finite weights, determinant, shape test, inverse or new MC/DP."
    )
    result = {"schema": "matching-one/p334-mixed-source-curvature-joint/v1",
              "allocation_commit": "fef79403c1727f7ba04273ac52c27e1e6696c248", "previous_shared_commit": OLD,
              "hessian_commit": HESSIAN, "rectangle_commit": RECTANGLE, "first_response_commit": NEW,
              "source_sha256": hashes, "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "p_ref": hs["p_ref"], "delta_cos4": DELTA, "sizes": sizes,
              "new_samples": 0, "new_DP": 0, "raw_fork_reads": 0, "finite_weight_evaluations": 0,
              "factor_alignment": "Hessian.LOO_factor=-Hessian.factor; rectangle.factor and previous.factor already LOO-centered",
              "boundary": boundary}
    (OUT / "score.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    lines += ["## Boundary", "", boundary, ""]
    (OUT / "REPORT.md").write_text("\n".join(lines))
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
