#!/usr/bin/env python3
"""Evaluate the new area-scale contrasts from already published Q velocities."""
from fractions import Fraction
import argparse
import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("do not replace an existing prediction packet")
    mp.mp.dps = 70
    result = {
        "schema": "matching-one.closed-source-weak-colour-predictions.v1",
        "source_commit": "85fd492312b597b3fa102ea913e4bcc7aeae2acf",
        "old_input": "notes/q-velocity-spin4-spectroscopy.md; velocities imported, oracle not rerun",
        "new_chain_rule": "Q_colour=exp(2t); area derivative exponent -x/2; target=-dx/dQ at Q=1",
        "statistic": "[V_(cN)/U_(cN)-V_N/U_N]/log(c), t=0, c is area ratio",
        "conditional_scope": "same critical Potts-Q continuation, nonzero single-field thermal-slope overlap, controlled corrections and projector derivative retained",
        "models": {},
    }
    for name, dimension, numerator in (
        ("thermal_Q4_epsilon", Fraction(21, 4), 9),
        ("four_leg_V22", Fraction(17, 4), 5),
    ):
        exponent = Fraction(13, 8)-(dimension-2)/2
        target = numerator*mp.sqrt(3)/(16*mp.pi)
        result["models"][name] = {
            "dimension_at_Q1": str(dimension), "U_N_exponent_at_Q1": str(exponent),
            "target_symbolic": f"{numerator}*sqrt(3)/(16*pi)",
            "target_decimal": mp.nstr(target, 60),
            "area_dilations": {
                str(c): {
                    "delta_R_symbolic": f"{numerator}*sqrt(3)*log({c})/(16*pi)",
                    "delta_R_decimal": mp.nstr(target*mp.log(c), 60),
                    "baseline_U_ratio_limit": mp.nstr(mp.mpf(c)**(mp.mpf(exponent.numerator)/exponent.denominator), 60),
                } for c in (2, 4)
            },
        }
    result.update({
        "target_gap_symbolic": "sqrt(3)/(4*pi)",
        "target_gap_decimal": mp.nstr(mp.sqrt(3)/(4*mp.pi), 60),
        "target_ratio": "9/5",
        "Q4_branch_boundary_t": "log(2)",
        "new_samples": 0, "new_finite_coupling_scores": 0, "velocity_oracle_calls": 0,
        "python": sys.executable,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"targets": {name: row["target_decimal"] for name, row in result["models"].items()},
                      "gap": result["target_gap_decimal"]}, indent=2))


if __name__ == "__main__":
    main()
