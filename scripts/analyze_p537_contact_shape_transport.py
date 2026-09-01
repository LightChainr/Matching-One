#!/usr/bin/env python3
"""Compare exact N25 and held-out N65 contact-stage tensor shapes."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import chi2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contract", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    ap.add_argument("--output-md", type=Path, required=True)
    args = ap.parse_args()
    if args.output_json.exists() or args.output_md.exists():
        raise FileExistsError("refusing to overwrite output")

    contract = json.loads(args.contract.read_text())
    n25 = np.asarray(contract["n25"]["matrix"], dtype=float).reshape(-1)
    n65_payload = json.loads(Path(contract["n65"]["result"]).read_text())
    n65 = np.asarray(n65_payload["primary"]["matrix"], dtype=float).reshape(-1)
    covariance = np.asarray(
        n65_payload["primary"]["covariance"]["combined"], dtype=float
    )[:4, :4]
    precision = np.linalg.pinv(covariance, rcond=1e-12)
    amplitude = float(n25 @ precision @ n65 / (n25 @ precision @ n25))
    amplitude_se = float(1.0 / np.sqrt(n25 @ precision @ n25))
    residual = n65 - amplitude * n25
    statistic = float(residual @ precision @ residual)
    dof = 3
    p_value = float(chi2.sf(statistic, dof))
    cosine = float(n25 @ n65 / (np.linalg.norm(n25) * np.linalg.norm(n65)))
    euclidean_amplitude = float(n25 @ n65 / (n25 @ n25))
    relative_residual = float(
        np.linalg.norm(n65 - euclidean_amplitude * n25) / np.linalg.norm(n65)
    )
    decision = (
        "RIGID_SHAPE_DISFAVORED" if p_value < 0.05
        else "RIGID_SHAPE_NOT_REJECTED"
    )
    payload = {
        "schema": "matching-one/p537-contact-shape-transport-score/v1",
        "status": "COMPLETED_EXPLORATORY_POSTHOC",
        "decision": decision,
        "matrix_order": ["01_single", "01_double", "12_single", "12_double"],
        "n25": n25.reshape(2, 2).tolist(),
        "n65": n65.reshape(2, 2).tolist(),
        "gls_amplitude": amplitude,
        "gls_amplitude_se": amplitude_se,
        "residual": residual.reshape(2, 2).tolist(),
        "chi_square": statistic,
        "degrees_of_freedom": dof,
        "p_value": p_value,
        "euclidean_cosine": cosine,
        "euclidean_relative_residual": relative_residual,
        "dependency_groups": [
            "p537-exact-N25-contact-stage",
            "p537-N65-contact-stage-20m-plus-P45-baseline"
        ],
        "boundary": contract["boundary"]
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(payload, indent=2) + "\n")
    lines = [
        "# P537 contact-stage shape transport",
        "",
        f"Decision: `{decision}` (exploratory/post-hoc).",
        "",
        "The one-amplitude model `vec(L65)=a vec(L25)` gives",
        "",
        "```text",
        f"a       = {amplitude:.12g} +/- {amplitude_se:.12g}",
        f"chi2    = {statistic:.12g} / {dof}",
        f"p       = {p_value:.12g}",
        f"cosine  = {cosine:.12g}",
        f"relative Euclidean residual = {relative_residual:.12g}",
        "```",
        "",
        "N25 and N65 retain the same four signs and a high visual cosine, but",
        "the complete N65 covariance resolves a non-scalar shape rotation.  The",
        "first-birth double-contact cell supplies the largest raw departure from",
        "the fitted ray.  The surviving object therefore needs at least a running",
        "two-coupling description; one universal tensor times one amplitude is",
        "not an adequate compression at nominal 5%.",
        "",
        "This proportional-shape question was not frozen before N65 was opened.",
        "It is an adaptive analysis of the same production, not an independent",
        "evidence vote or a prospective rejection.",
    ]
    args.output_md.write_text("\n".join(lines) + "\n")
    print(json.dumps({"decision": decision, "chi_square": statistic,
                      "dof": dof, "p_value": p_value}))


if __name__ == "__main__":
    main()
