#!/usr/bin/env python3
"""Exact decimation transport of published UNMARKED U values and covariance.

No configurations, source marks, fresh validation data, or fitted coefficients.
The numerical statements are old-baseline endpoint predictions, not new MC.
"""
import argparse
import hashlib
import json
import subprocess
from fractions import Fraction
from pathlib import Path

import numpy as np

REV = "7da1eeb0e51cf430987dbf204d23713c2ab5a46c"
SOURCE = "results/norm4-source-endpoint-1m/latest.json"
NS = (65, 85, 130, 170, 260, 340)
PAIRS = ((130, 65), (170, 85), (260, 130), (340, 170))


def h4(z):
    a, b = z
    return Fraction(a**4 - 6*a*a*b*b + b**4, (a*a+b*b)**2)


def canonical(z):
    return tuple(sorted(map(abs, z), reverse=True))


def period_map(parent, child):
    mapped = {}
    canonical_child = {canonical(child[g]): g for g in ("first", "second")}
    for g in ("first", "second"):
        a, b = parent[g]
        assert (a+b) % 2 == 0 and (b-a) % 2 == 0
        direct = ((a+b)//2, (b-a)//2)
        target = canonical_child[canonical(direct)]
        assert h4(direct) == -h4((a, b)) == h4(child[target])
        mapped[g] = {"parent_gaussian": [a, b], "divide_by_1_plus_i": direct,
                     "canonical_child_orientation": target,
                     "parent_cos4": str(h4((a, b))), "child_cos4": str(h4(direct))}
    assert {r["canonical_child_orientation"] for r in mapped.values()} == {"first", "second"}
    return mapped


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    raw = subprocess.check_output(["git", "show", f"{REV}:{SOURCE}"])
    src = json.loads(raw)
    idx = [src["labels"].index(f"N{n}.U_anchor") for n in NS]
    u = np.array([src["estimates"][f"N{n}.U_anchor"]["value"] for n in NS])
    cov = np.asarray(src["covariance"])[np.ix_(idx, idx)]
    factor = 2.0**(13/8)
    # A single deterministic map: six native U, four saturated U, four differences.
    transform = np.zeros((14, 6))
    transform[:6] = np.eye(6)
    labels = [f"N{n}.native_U" for n in NS]
    labels += [f"N{n}.saturated_U" for n, _ in PAIRS]
    labels += [f"N{n}.saturated_minus_native" for n, _ in PAIRS]
    for i, (parent, child) in enumerate(PAIRS):
        transform[6+i, NS.index(child)] = factor
        transform[10+i] = transform[6+i]
        transform[10+i, NS.index(parent)] -= 1
    values = transform @ u
    joint_cov = transform @ cov @ transform.T
    groups, inherited_cov = {}, np.zeros((6, 6))
    for name, group in src["covariance_contributions"].items():
        selected = np.asarray(group["covariance"])[np.ix_(idx, idx)]
        inherited_cov += selected
        if group["stage"] == "marked_source":
            assert np.max(np.abs(selected)) < 1e-25
            continue  # these marks contribute no uncertainty to unmarked anchors.
        loo = np.asarray(group["delete_one_vectors"])[:, idx]
        groups[name] = {"Ns": group["Ns"], "batch_counts": group["batch_counts"],
                        "delete_one_batch_ids": group["delete_one_batch_ids"],
                        "operation": group["operation"],
                        "selected_native_covariance": selected.tolist(),
                        "selected_native_delete_one_vectors": loo.tolist(),
                        "joint_covariance": (transform @ selected @ transform.T).tolist()}
    assert np.allclose(inherited_cov, cov, rtol=1e-12, atol=1e-15)
    estimates = {label: {"value": float(value), "se": float(se)}
                 for label, value, se in zip(labels, values, np.sqrt(np.diag(joint_cov)))}
    rows = []
    for i, (parent, child) in enumerate(PAIRS):
        p_root = src["by_N"][str(child)]["anchor"]["root"]
        rows.append({"N_parent": parent, "N_child": child,
                     "period_map": period_map(src["unmarked_inputs"][str(parent)]["design"],
                                               src["unmarked_inputs"][str(child)]["design"]),
                     "p_white_endpoint": 1-p_root,
                     "p_black_endpoint": 1,
                     "native": estimates[labels[NS.index(parent)]],
                     "saturated": estimates[labels[6+i]],
                     "difference": estimates[labels[10+i]]})
    output = {"status": "EXACT_ENDPOINT_MAP_WITH_INHERITED_OLD_BASELINE_UNCERTAINTY",
              "family": "p_A=s+(1-s)*p, p_B=p; thermal derivative holds s fixed",
              "identity": "U_parent(s=1)=2^(13/8)*U_child(s=0)",
              "source_commit": REV, "source_path": SOURCE,
              "source_sha256": hashlib.sha256(raw).hexdigest(),
              "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "selected_input_labels": [src["labels"][i] for i in idx],
              "unmarked_input_provenance": src["unmarked_inputs"],
              "factor_exact": "2^(13/8)", "factor_numeric": factor,
              "labels": labels, "estimates": estimates,
              "linear_transform": transform.tolist(), "joint_covariance": joint_cov.tolist(),
              "dependency_groups": groups, "rows": rows,
              "new_samples": 0, "configuration_replays": 0,
              "fresh_validation_data_used": False, "source_coefficients_fitted": False,
              "boundary": "Exact finite-topology transport; old unmarked numerical anchors and their existing covariance. No fresh stochastic endpoint measurement, intermediate-s curve, curvature sign, continuum field identity, or new independent evidence is claimed."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(output, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
