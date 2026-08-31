#!/usr/bin/env python3
"""Final shared-LOO handoff: lifetime, full A/E cells, shape energy and source debt."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUTS = {
    "lifetime": ("be31a113", "results/p334-lifetime-square-mechanism/score.json"),
    "shape": ("da0080ec", "results/p334-centered-birth-shape/batch_vectors.json"),
    "nine_cells": ("bb79fd47", "results/p334-nine-layer-complete-ae/batch_vectors.json"),
    "connected": ("9ed1e5082ac114d818da03d07e6cf2a315d75023", "results/p334-first-birth-winner-connected-coupling/score.json"),
}
GROUPS = {"00": [0], "11": [4], "22": [8], "01+10": [1, 3], "02+20": [2, 6], "12+21": [5, 7]}
OUT = ROOT/"results/p334-complete-shape-source-joint"


def linear_loo(batches):
    return (20*batches.mean(axis=0)-batches)/19


def covariance(loo):
    centered = loo-loo.mean(axis=0)
    return 19/20*centered.T@centered


def shape_coordinates(raw, lifetime, life_labels, n, delta):
    center, width, canonical = (raw[0]-raw[4])/delta, (raw[1]-raw[5])/delta, (raw[2]-raw[6])/delta
    ensemble = lifetime[life_labels.index("H4.mixture_Y_variance")]
    center_spread = lifetime[life_labels.index("H4.C_variance")]/(n+1)**2
    connected = lifetime[life_labels.index("H4.K1K2_covariance")]/(n+1)**2
    return np.array([center, width, center+width, canonical, ensemble,
                     center+width-ensemble, center-width, connected, center_spread])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    data, provenance = {}, {}
    for name, (commit, path) in INPUTS.items():
        commit = subprocess.check_output(["git", "rev-parse", commit+"^{commit}"], cwd=ROOT, text=True).strip()
        blob = subprocess.check_output(["git", "show", commit+":"+path], cwd=ROOT)
        data[name] = json.loads(blob)
        provenance[name] = {"commit": commit, "path": path, "sha256": sha256(blob).hexdigest()}
    source = data["lifetime"]["source_commit"]
    if (source != data["shape"]["source_commit"] or source != data["nine_cells"]["full_birth_commit"]
            or source != data["connected"]["birth_commit"]):
        raise ValueError("Different original full-birth dependency blocks")
    shape_names = ["rank_center_energy", "rank_lifetime_energy", "R_ref_total", "Q_ref_canonical",
                   "ensemble_Y_variance", "mixture_mean_bias_squared", "joint_center_minus_lifetime",
                   "connected_K1K2_normalized", "ensemble_center_variance"]
    sizes = {}
    for n in (325, 425):
        w, s, a, c = (data[k]["sizes"][str(n)] for k in INPUTS)
        if any(r["batch_ids"] != list(range(20)) for r in (w, s, a, c)):
            raise ValueError("Original paired batch IDs differ")
        if (w["batch_denominators"] != [1000]*20 or s["batch_denominators"] != [1000]*20
                or a["samples_per_batch"] != 1000 or c["batch_counter_count"] != 1000):
            raise ValueError("Original complete batch denominators differ")
        delta = w["delta_cos4"]
        if a["delta_cos4"] != delta or c["delta_cos4"] != delta:
            raise ValueError("H4 conventions differ")
        path = f"results/p334-full-birth-archive/N{n}.csv"
        if s["source_sha256"] != c["birth_sha256"] or s["source_sha256"] != data["nine_cells"]["source_sha256"][path]:
            raise ValueError("Original CSV hashes differ")
        wb, sb = np.array(w["joint_LOO_vectors"]), np.array(s["joint_20_batch_means"])
        ab = np.array(a["joint_20_batch_means"])
        cb = np.array(c["leave_one_original_batch_out"])
        sl, al = linear_loo(sb), linear_loo(ab)
        wp, sp = np.array(w["estimate"]), sb.mean(axis=0)
        ic, iw = (w["labels"].index("H4."+name) for name in ("C_mean", "W_mean"))
        birth_mean = np.array([wp[ic]-wp[iw]/2, wp[ic]+wp[iw]/2])
        birth_mean_loo = np.column_stack((wb[:, ic]-wb[:, iw]/2, wb[:, ic]+wb[:, iw]/2))
        sh = shape_coordinates(sp, wp, w["labels"], n, delta)
        shl = np.array([shape_coordinates(sl[b], wb[b], w["labels"], n, delta) for b in range(20)])
        total = ab.reshape(20, 9, 9)[:, :, 1:].sum(axis=1)
        if np.max(np.abs(total-np.array(a["full_AE_total_batch_means"]))) > 1e-12:
            raise ValueError("Nine cells do not add to full A/E")
        group_blocks = [ab.reshape(20, 9, 9)[:, indices, :].sum(axis=1) for indices in GROUPS.values()]
        group_b = np.column_stack(group_blocks)
        group_labels = ["group."+name+"."+field for name in GROUPS for field in data["nine_cells"]["cell_fields"]]
        labels = (["lifetime."+name for name in w["labels"]]
                  +["birth_mean.H4.K1", "birth_mean.H4.K2"]
                  +["raw_shape."+name for name in data["shape"]["labels"]]
                  +["shape.H4."+name for name in shape_names]
                  +["global."+name for name in data["nine_cells"]["effect_fields"]]
                  +group_labels+["source."+name for name in data["connected"]["derived_labels"]])
        point = np.r_[wp, birth_mean, sp, sh, total.mean(axis=0), group_b.mean(axis=0), c["derived_point"]]
        joint = np.column_stack((wb, birth_mean_loo, sl, shl, linear_loo(total), linear_loo(group_b), cb))
        cov = covariance(joint)
        cell_mean, cell_cov = ab.mean(axis=0), covariance(al)
        raw_second = np.zeros((81, 81))
        for k, block in enumerate(a["cell_individual_raw_second_moments"]):
            raw_second[k*9:(k+1)*9, k*9:(k+1)*9] = block
        individual_cov = raw_second-np.outer(cell_mean, cell_mean)
        # Distinct joint cells are disjoint. Raw cross moments vanish, but their
        # covariance is -outer(means), not zero and not an independence license.
        g_indices = [labels.index("global."+field) for field in data["nine_cells"]["effect_fields"]]
        shape_indices = [labels.index("shape.H4."+name) for name in shape_names]
        shape_cov = cov[np.ix_(shape_indices, shape_indices)]
        closures = {
            "global_baseline_E_integral_vs_W_mean": float(point[labels.index("global.baseline.p_integral.E")]+wp[w["labels"].index("H4.W_mean")]/(n+1)),
            "Rref_equals_ensemble_plus_mean_bias": float(sh[2]-sh[4]-sh[5]),
            "rank_lifetime_vs_raw_W2": float(sh[1]-wp[w["labels"].index("H4.W_second_moment")]/(4*(n+1)**2)),
            "connected_D_plus_G": float(c["derived_point"][data["connected"]["derived_labels"].index("H4_direct_connected_debt")]
                                         +c["derived_point"][data["connected"]["derived_labels"].index("H4_collective_connected_debt")]),
            "shared_F2_shift_pref_max_batch": float(np.max(np.abs((total[:, 4]-total[:, 0])-(total[:, 5]-total[:, 1])))),
            "shared_F2_shift_integral_max_batch": float(np.max(np.abs((total[:, 6]-total[:, 2])-(total[:, 7]-total[:, 3])))),
        }
        sizes[str(n)] = {"delta_cos4": delta, "batch_ids": list(range(20)), "samples_per_batch": 1000,
            "labels": labels, "estimate": point.tolist(), "se": np.sqrt(np.diag(cov)).tolist(),
            "joint_LOO_vectors": joint.tolist(), "joint_covariance": cov.tolist(), "joint_rank_at_most": 19,
            "identities_residual": closures,
            "global_AE": {"labels": data["nine_cells"]["effect_fields"], "mean": point[g_indices].tolist(),
                          "se": np.sqrt(np.diag(cov)[g_indices]).tolist(), "covariance": cov[np.ix_(g_indices, g_indices)].tolist(),
                          "joint_20_batch_means": total.tolist()},
            "shape_H4": {"labels": shape_names, "mean": sh.tolist(), "se": np.sqrt(np.diag(shape_cov)).tolist(), "covariance": shape_cov.tolist()},
            "six_transpose_groups": {"group_cells": GROUPS, "labels": group_labels, "joint_20_batch_means": group_b.tolist()},
            "nine_cells": {"labels": data["nine_cells"]["labels"], "joint_20_batch_means": ab.tolist(),
                           "mean": cell_mean.tolist(), "mean_se": np.sqrt(np.diag(cell_cov)).tolist(),
                           "mean_covariance": cell_cov.tolist(), "individual_covariance": individual_cov.tolist(),
                           "individual_covariance_convention": "ddof0; block-diagonal raw second moments minus full outer mean; cross-cell covariances retained"}}
    result = {"schema": "matching-one/p334-complete-shape-source-joint/v1", "sources": provenance,
        "full_birth_commit": source, "sizes": sizes, "new_MC": 0, "new_DP": 0, "new_path_replay": 0,
        "boundary": "All coordinates reuse the same paired block per N. Nonlinear lifetime/source coordinates use the supplied population-recentered LOO; other coordinates use the same original deleted batch. Nine joint-rank cells are not the four R1-flag states. Full covariance is preserved without inversion or a new omnibus test."}
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Complete A/E, thermal shape and source-connected debt", ""]
    for n, r in sizes.items():
        lines += [f"## N{n}", "", "| Global/shape readout | Mean | Shared-batch SE |", "|---|---:|---:|"]
        for name, value, se in zip(r["labels"], r["estimate"], r["se"]):
            if name.startswith(("global.", "shape.H4.", "source.H4_", "birth_mean.H4.")):
                lines.append(f"| {name} | {value:.10g} | {se:.6g} |")
        lines += ["", "| Transpose group | safe pref A | safe pref E | safe integral A | safe integral E |", "|---|---:|---:|---:|---:|"]
        for group in GROUPS:
            entries = []
            for ep in ("p_ref", "p_integral"):
                for observer in ("A", "E"):
                    ix = r["labels"].index(f"group.{group}.safe.{ep}.{observer}")
                    entries.append(f"{r['estimate'][ix]:.8g} +/- {r['se'][ix]:.5g}")
            lines.append("| "+group+" | "+" | ".join(entries)+" |")
        lines.append("")
    lines += [result["boundary"], ""]
    (args.output/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
