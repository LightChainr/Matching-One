#!/usr/bin/env python3
"""Score the frozen N25-to-N65 contact-stage scaling fingerprint."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import chi2


def collapse_n25(data):
    masks = data["contact_decomposition"]["matrices"]
    signed = np.zeros((2, 2))
    exposure = np.zeros((2, 2))
    for stage in range(2):
        for mask in ("1", "2"):
            signed[stage, 0] += sum(
                cell["midpoint"] for cell in masks[mask]["cells"][stage]
            )
            exposure[stage, 0] += sum(
                cell["midpoint"] for cell in masks[mask]["positive_mass"][stage]
            )
        signed[stage, 1] = sum(
            cell["midpoint"] for cell in masks["3"]["cells"][stage]
        )
        exposure[stage, 1] = sum(
            cell["midpoint"] for cell in masks["3"]["positive_mass"][stage]
        )
    return signed, exposure


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n25", required=True, type=Path)
    ap.add_argument("--n65", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)

    n25, e25 = collapse_n25(json.loads(args.n25.read_text()))
    d65 = json.loads(args.n65.read_text())
    n65 = np.asarray(d65["primary"]["matrix"], dtype=float)
    cov = np.asarray(d65["primary"]["matrix_covariance"], dtype=float)
    e = np.asarray(d65["positive_exposure"]["estimate"], dtype=float)
    e65 = np.column_stack((e[:, 0] + e[:, 1], e[:, 2]))
    ratio = 65 / 25
    powers = np.array([3, 29 / 8, 3, 3], dtype=float)
    predicted = n25.ravel() * ratio ** (-powers)

    def score(prediction, df):
        residual = n65.ravel() - prediction
        q = float(residual @ np.linalg.solve(cov, residual))
        return {
            "prediction": prediction.reshape(2, 2).tolist(),
            "Q": q,
            "df": df,
            "p_value": float(chi2.sf(q, df)),
            "standardized_marginal_residual": (
                residual / np.sqrt(np.diag(cov))
            ).tolist(),
        }

    split = score(predicted, 4)
    common_n3 = score(n25.ravel() * ratio ** -3, 4)
    inv = np.linalg.inv(cov)
    anchor = n25.ravel()
    lam = float(anchor @ inv @ n65.ravel() / (anchor @ inv @ anchor))
    free_amplitude = score(lam * anchor, 3)
    free_amplitude["lambda"] = lam

    chi25 = float(n25[0, 0] * n25[1, 1] / (n25[0, 1] * n25[1, 0]))
    chi65 = float(d65["primary"]["projective_cross_ratio"])
    chi65_se = float(d65["primary"]["projective_cross_ratio_se"])
    invariant25 = -chi25 * 25 ** (-5 / 8)
    invariant65 = -chi65 * 65 ** (-5 / 8)
    invariant65_se = chi65_se * 65 ** (-5 / 8)
    density25 = n25 / e25
    density65 = n65 / e65

    def exponent(x25, x65, sign=-1):
        return sign * math.log(abs(x65 / x25)) / math.log(ratio)

    payload = {
        "schema": "matching-one/p537-contact-stage-scale/v1",
        "status": "POST_N65_TWO_SCALE_MECHANISM_FINGERPRINT",
        "row_order": ["01", "12"],
        "column_order": ["single", "double"],
        "N25_matrix": n25.tolist(),
        "N65_matrix": n65.tolist(),
        "split_power_model": {
            "powers_in_N_row_major": powers.tolist(),
            **split,
        },
        "common_N_minus_3": common_n3,
        "free_common_amplitude": free_amplitude,
        "projective_invariant": {
            "definition": "N^(-5/8)*(-K01s*K12d/(K01d*K12s))",
            "N25": invariant25,
            "N65": invariant65,
            "N65_se": invariant65_se,
            "relative_center_difference": invariant65 / invariant25 - 1,
            "chi25": chi25,
            "chi65": chi65,
            "chi65_se": chi65_se,
        },
        "entry_double_decomposition": {
            "N25_exposure": e25[0, 1],
            "N65_exposure": e65[0, 1],
            "exposure_decay_power": exponent(e25[0, 1], e65[0, 1]),
            "N25_conditional_density": density25[0, 1],
            "N65_conditional_density": density65[0, 1],
            "conditional_density_growth_power": exponent(
                density25[0, 1], density65[0, 1], sign=1
            ),
            "signed_mass_decay_power": exponent(n25[0, 1], n65[0, 1]),
        },
        "interpretation": "three N^-3 cells plus one entry-double cell with an additional N^-5/8 thermal suppression in conditional signed strength",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({"split": split, "common_N3": common_n3}, indent=2))


if __name__ == "__main__":
    main()
