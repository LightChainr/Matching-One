#!/usr/bin/env python3
"""Read a fixed-support lower-bound parameter with its sampling uncertainty."""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import numpy as np
from p437_fixed_support_mc import DENOMINATOR, energy_numerator
from p437_positive_difference_bridge import H5


def score(directory):
    metadata = json.loads((directory / "run.json").read_text())
    path = directory / "batches.json"
    if hashlib.sha256(path.read_bytes()).hexdigest() != metadata["batch_sha256"]:
        raise ValueError("batch hash mismatch")
    rows = json.loads(path.read_text())
    classes = Counter()
    vectors = []
    for row in rows:
        vector = np.zeros(8)
        for value in row["classes"]:
            a, b, c = value["child_difference_numerators"]
            count = value["count"]
            numerator = energy_numerator((a, b, c))
            if numerator != value["energy_numerator"]:
                raise AssertionError("integer energy reconstruction differs")
            classes[(a, b, c)] += count
            vector += count * np.array([numerator / DENOMINATOR, a*a/1024, b*b/1024, c*c/1024,
                                        a*b/1024, a*c/1024, b*c/1024, int(numerator != 0)])
        vectors.append(vector / row["samples"])
    x = np.array(vectors)
    mean = x.mean(axis=0)
    cov = np.cov(x, rowvar=False, ddof=1) / len(x)
    samples = sum(row["samples"] for row in rows)
    nonzero = sum(count for value, count in classes.items() if energy_numerator(value))
    total_numerator = sum(count * energy_numerator(value) for value, count in classes.items())
    energy = Fraction(total_numerator, DENOMINATOR * samples)
    phase_classes = []
    for value, count in sorted(classes.items()):
        phase_classes.append({"child_difference_numerators": value, "count": count,
                              "squared_derivative_exact": str(Fraction(energy_numerator(value), DENOMINATOR))})
    return {"schema": "matching-one/p437-fixed-support-score/v1", "samples": samples,
            "B_S": {"mean_exact": str(energy), "mean": float(energy), "batch_se": float(np.sqrt(cov[0,0]))},
            "h5_B_S_population_lower_bound_parameter": {"estimate_exact": str(H5 * energy),
                "estimate": float(H5 * energy), "batch_se": float(H5) * float(np.sqrt(cov[0,0])),
                "not_a_statistically_certain_lower_bound": True},
            "nonzero": {"count": nonzero, "rate": nonzero/samples,"batch_se": float(np.sqrt(cov[7,7])),
                "batches_with_nonzero": sum(row["nonzero"] > 0 for row in rows)},
            "coordinate_order": ["B_S","D0_energy","D1_energy","D2_energy","D0D1","D0D2","D1D2","nonzero_rate"],
            "point": mean.tolist(), "full_covariance_of_mean": cov.tolist(), "value_classes": phase_classes,
            "runtime": {key: metadata[key] for key in ("wall_seconds","cpu_seconds")},
            "decision": "fixed_support_readable" if nonzero else "no_nonzero_in_fixed20k_not_population_zero",
            "boundary": "Fixed support chosen by prior exact certificate, no scan. Localized positive energy endpoint, not unbiased total high-pass reconstruction. Mean and SE do not turn the population inequality into a certain numerical bound."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(score(args.directory), indent=2, sort_keys=True, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
