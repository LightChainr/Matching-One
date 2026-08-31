#!/usr/bin/env python3
"""Original-batch continuous-clock dipole and connected-response decomposition."""
import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OLD_COMMIT = "b582015e64e2d8a59e591c4822b14dedaea58b0f"
OLD_DIR = "results/p334-common-label-tangent-joint"
THERMAL_COMMIT = "9059776d866287e0cdb95e0a5f079843905cbeb9"
THERMAL_PATH = "results/p334-euler-thermal-moments/score.json"
MOMENT_COMMIT = "f4682eb379b5709a2840faf92beef44ff27f6f23"
MOMENT_PATH = "results/p334-continuous-center-lifetime-moments/batch_moments.json"
PLATEAU_COMMIT = "1e8549b5b43c27fc1f2d6691caed23bf75cf0d15"
PLATEAU_PATH = "results/p334-plateau-shape-tangent/score.json"
OUT = ROOT/"results/p334-euler-dipole-connected-clock"
ORIENTATIONS = ("first", "second")
MARKS = ("plus", "minus")
MOMENTS = ("C", "W", "CW", "C2", "W2")


def derive(baseline, tangent, n, delta, pref):
    """Baseline shape (orientation,5), H shape (mark,orientation,5)."""
    values = {}
    orientation_values = []
    for oi, orientation in enumerate(ORIENTATIONS):
        mc, mw, mcw, mc2, mw2 = baseline[oi]
        vc, vw, ccw = mc2-mc*mc, mw2-mw*mw, mcw-mc*mw
        v1, v2, c12 = vc+vw/4-ccw, vc+vw/4+ccw, vc-vw/4
        m12, mu1 = mc2-mw2/4, mc-mw/2
        orderstat_cov = (mu1-m12)/(n+1)
        if v1 <= 0 or v2 <= 0:
            raise ValueError("Continuous correlation denominator is not identified")
        root_v = np.sqrt(v1*v2); rho = c12/root_v
        base = {"mu_C": mc, "mu_W": mw, "var_C": vc, "var_W": vw, "cov_CW": ccw,
                "var_tau1": v1, "var_tau2": v2, "cov_tau12": c12, "rho_tau12": rho,
                "orderstat_cov_tau12": orderstat_cov, "intrinsic_rank_cov_tau12": c12-orderstat_cov}
        values.update({f"baseline.{orientation}.{k}": v for k, v in base.items()})
        mark_values = []
        for mi, mark in enumerate(MARKS):
            hc, hw, hcw, hc2, hw2 = tangent[mi, oi]
            dvc, dvw = hc2-2*mc*hc, hw2-2*mw*hw
            dccw = hcw-mc*hw-mw*hc
            dv1, dv2, dc12 = dvc+dvw/4-dccw, dvc+dvw/4+dccw, dvc-dvw/4
            drho = dc12/root_v-rho*(dv1/v1+dv2/v2)/2
            orderstat = ((hc-hw/2)-(hc2-hw2/4))/(n+1)
            row = {"delta_var_C": dvc, "delta_var_W": dvw, "delta_cov_CW": dccw,
                   "delta_var_tau1": dv1, "delta_var_tau2": dv2,
                   "delta_cov_tau12": dc12, "delta_rho_tau12": drho,
                   "delta_orderstat_cov_tau12": orderstat, "delta_intrinsic_rank_cov_tau12": dc12-orderstat,
                   "dipole_center_displacement": -mw*hc,
                   "dipole_mean_lifetime": -(mc-pref)*hw,
                   "dipole_connected_spread": -dccw,
                   "dipole_total": -hcw+pref*hw}
            values.update({f"{mark}->{orientation}.{k}": v for k, v in row.items()})
            mark_values.append(row)
        orientation_values.append(mark_values)
    # Only now transform the derived orientation quantities into S and D.
    for mi, mark in enumerate(MARKS):
        first, second = orientation_values[0][mi], orientation_values[1][mi]
        for key in first:
            values[f"{mark}->S.{key}"] = (first[key]+second[key])/2
            values[f"{mark}->D.{key}"] = (first[key]-second[key])/delta
    return values


def main():
    hashes = {}
    def read(commit, path):
        blob = subprocess.check_output(["git", "show", commit+":"+path], cwd=ROOT)
        hashes[commit+":"+path] = sha256(blob).hexdigest()
        return blob
    old = json.loads(read(OLD_COMMIT, OLD_DIR+"/score.json"))
    thermal = json.loads(read(THERMAL_COMMIT, THERMAL_PATH))
    moment = json.loads(read(MOMENT_COMMIT, MOMENT_PATH))
    plateau = json.loads(read(PLATEAU_COMMIT, PLATEAU_PATH))
    if moment["source_commit"] != old["source_commit"] or moment["policy_commit"] != thermal["source_commit"]:
        raise ValueError("Different fork source or common-label policy")
    if plateau["source_commit"] != MOMENT_COMMIT:
        raise ValueError("Plateau readout is from a different continuous-moment source")
    OUT.mkdir(parents=True, exist_ok=True)
    sizes = {}
    for ns in ("325", "425"):
        n = int(ns); o, t, m, pl = old["sizes"][ns], thermal["sizes"][ns], moment["sizes"][ns], plateau["sizes"][ns]
        if o["batch_ids"] != m["batch_ids"] or t["batch_ids"] != m["batch_ids"] or pl["batch_ids"] != m["batch_ids"] or m["batch_ids"] != list(range(20)):
            raise ValueError("Original twenty-batch alignment differs")
        raw = np.array(m["joint_20_batch_means"]); mean = raw.mean(axis=0)
        raw_loo = (20*mean-raw)/19
        baseline_ix = [m["labels"].index(f"all.baseline.{ori}.{field}") for ori in ORIENTATIONS for field in MOMENTS]
        tangent_ix = [m["labels"].index(f"all.H.{mark}.{ori}.{field}") for mark in MARKS for ori in ORIENTATIONS for field in MOMENTS]
        def from_raw(row):
            return derive(row[baseline_ix].reshape(2, 5), row[tangent_ix].reshape(2, 2, 5), n, m["delta_cos4"], thermal["p_ref"])
        point = from_raw(mean); derived_labels = list(point)
        derived_loo = np.array([list(from_raw(row).values()) for row in raw_loo])
        thermal_raw = np.array(t["joint_20_batch_means"])
        thermal_mean = thermal_raw.mean(axis=0)
        thermal_loo = (20*thermal_mean-thermal_raw)/19
        new_loo = np.column_stack((raw_loo, thermal_loo, derived_loo))
        new_factor = np.sqrt(19/20)*(new_loo-new_loo.mean(axis=0))
        saved = json.loads(gzip.decompress(read(OLD_COMMIT, OLD_DIR+"/"+o["complete_covariance_factor_file"])))
        new_labels = (["continuous_moment."+k for k in m["labels"]]+["thermal_moment."+k for k in t["labels"]]
                      +["connected_clock."+k for k in derived_labels])
        labels = saved["labels"]+new_labels+["plateau."+k for k in pl["labels"]]
        factor = np.column_stack((np.array(saved["factor"]), new_factor, np.array(pl["factor"])))
        values = {"continuous_moment."+k: v for k, v in zip(m["labels"], mean)}
        values.update({"thermal_moment."+k: v for k, v in zip(t["labels"], thermal_mean)})
        values.update({"connected_clock."+k: v for k, v in point.items()})
        values.update({"plateau."+k: v for k, v in zip(pl["labels"], pl["estimate"])})
        focused = ["connected_clock."+k for k in derived_labels if k.startswith("baseline.")]
        for mark, axis in (("plus", "S"), ("minus", "D")):
            focused += ["connected_clock."+k for k in derived_labels if k.startswith(f"{mark}->{axis}.")]
            focused += [f"thermal_moment.all.{mark}->{axis}.E.{k}" for k in ("I0", "I1", "I2", "dipole_at_p_ref", "quadrupole_at_p_ref")]
            focused += [f"plateau.{mark}->{axis}.{k}" for k in ("mass", "centroid", "variance", "width", "centroid_minus_unweighted_C")]
        ix = [labels.index(k) for k in focused]
        cov = factor[:, ix].T@factor[:, ix]
        filename = f"N{n}.complete_common_factor.json.gz"
        packed = {"labels": labels, "factor": factor.tolist(), "batch_ids": list(range(20)),
            "convention": "factor.T@factor; same original twenty deleted-batch rows; no inverse",
            "new_raw_moment_labels": m["labels"], "new_raw_moment_20_batch_means": raw.tolist(),
            "thermal_moment_labels": t["labels"], "thermal_moment_20_batch_means": thermal_raw.tolist(),
            "derived_labels": derived_labels, "derived_LOO": derived_loo.tolist(),
            "plateau_derived_labels": pl["labels"], "plateau_derived_LOO": pl["LOO"]}
        blob = gzip.compress((json.dumps(packed, separators=(",", ":"), allow_nan=False)+"\n").encode(), mtime=0)
        (OUT/filename).write_bytes(blob)
        derived_factor = new_factor[:, raw.shape[1]+thermal_raw.shape[1]:]
        dcov = derived_factor.T@derived_factor
        sizes[ns] = {"batch_ids": list(range(20)), "delta_cos4": m["delta_cos4"],
            "labels": focused, "estimate": [values[k] for k in focused], "se": np.sqrt(np.diag(cov)).tolist(),
            "focused_covariance": cov.tolist(), "derived_labels": derived_labels, "derived_estimate": list(point.values()),
            "derived_se": np.sqrt(np.diag(dcov)).tolist(), "derived_covariance": dcov.tolist(),
            "derived_LOO": derived_loo.tolist(), "rho_denominator_status": "positive variance in both orientations and all twenty LOO replicates",
            "complete_covariance_factor_file": filename, "complete_covariance_factor_sha256": sha256(blob).hexdigest(),
            "complete_coordinate_count": len(labels)}
    result = {"schema": "matching-one/p334-euler-dipole-connected-clock/v1", "source_commit": moment["source_commit"],
        "common_label_policy_commit": moment["policy_commit"], "previous_shared_commit": OLD_COMMIT,
        "moment_commit": MOMENT_COMMIT, "thermal_moment_commit": THERMAL_COMMIT,
        "plateau_readout_commit": PLATEAU_COMMIT, "source_sha256": hashes,
        "p_ref": thermal["p_ref"], "sizes": sizes, "new_MC": 0, "new_DP": 0, "new_raw_fork_reads": 0,
        "estimand": "Orientation-level full-population continuous-clock covariance derivative; products formed before S/D and separately inside every pooled leave-one-batch replicate",
        "identity": "dipole=-muW*HC-(muC-pref)*HW-deltaCovCW; deltaCovCW=(deltaVarTau2-deltaVarTau1)/2; deltaCovTau12=deltaVarC-deltaVarW/4",
        "orderstat_identity": "orderstat_cov=(muTau1-Etau1tau2)/(N+1); intrinsic_rank_cov=total_cov-orderstat_cov, all continuous scale",
        "boundary": "Thermal dipole does not identify a joint copula; CovCW is a marginal variance imbalance. CovTau12 additionally contains shared order-statistic timing, separated explicitly here. Full same-source covariance is retained without fitting, inverse covariance, new sampling or independent-evidence claims."}
    (OUT/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Thermal dipole and continuous center-lifetime response", ""]
    for ns, r in sizes.items():
        lines += [f"## N{ns}", "", "| Component | plus -> S | minus -> D |", "|---|---:|---:|"]
        for key in ("dipole_center_displacement", "dipole_mean_lifetime", "dipole_connected_spread", "dipole_total",
                    "delta_cov_CW", "delta_var_C", "delta_var_W", "delta_cov_tau12", "delta_orderstat_cov_tau12",
                    "delta_intrinsic_rank_cov_tau12", "delta_rho_tau12"):
            cells = []
            for mark, axis in (("plus", "S"), ("minus", "D")):
                i = r["derived_labels"].index(f"{mark}->{axis}.{key}")
                cells.append(f"{r['derived_estimate'][i]:.10g} +/- {r['derived_se'][i]:.6g}")
            lines.append("| "+key+" | "+" | ".join(cells)+" |")
        lines.append("")
        lines += ["| Supplied normalized plateau coordinate | plus -> S | minus -> D |", "|---|---:|---:|"]
        for key in ("centroid", "variance", "centroid_minus_unweighted_C"):
            cells = []
            for mark, axis in (("plus", "S"), ("minus", "D")):
                i = r["labels"].index(f"plateau.{mark}->{axis}.{key}")
                cells.append(f"{r['estimate'][i]:.10g} +/- {r['se'][i]:.6g}")
            lines.append("| "+key+" | "+" | ".join(cells)+" |")
        lines.append("")
    lines += [result["boundary"], ""]
    (OUT/"REPORT.md").write_text("\n".join(lines)); print("\n".join(lines))


if __name__ == "__main__": main()
