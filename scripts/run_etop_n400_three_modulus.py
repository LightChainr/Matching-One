#!/usr/bin/env python3
"""Acquire the new homothetic scale with the existing threshold-rank engine."""
import argparse
import json

import numpy as np

import run_etop_n100_three_modulus as engine

ROOT = engine.ROOT
OUT = ROOT/"results/etop-n400-three-modulus"
CONTRACT = ROOT/"experiments/etop_n400_three_modulus_20260831.json"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--freeze", required=True)
    args = parser.parse_args()
    contract = json.loads(CONTRACT.read_text())
    # This process runs only N400; reuse the old launcher/reader, not its N50
    # calibration or N100-specific report. No old production is touched.
    engine.OUT, engine.FREEZE = OUT, args.freeze
    if args.run:
        engine.launch(contract)
    batches = np.column_stack([engine.read_batch_vectors(contract, s) for s in contract["shapes"]])
    mean = batches.mean(axis=0)
    covariance = np.cov(batches, rowvar=False, ddof=1)/len(batches)
    result = {"contract": contract, "prediction_freeze_commit": args.freeze,
        "shape_order": [s["name"] for s in contract["shapes"]], "field_order": engine.FIELDS,
        "mean": mean.tolist(), "covariance": covariance.tolist(), "batch_vectors": batches.tolist(),
        "same_area_models": engine.same_area_models(mean, covariance),
        "raw_sha256": {str(p.relative_to(ROOT)): engine.digest(p) for p in sorted((OUT/"raw").iterdir())},
        "inference": "N400 is independent of N100; all N400 shapes are one shared-counter block. Fixed-p comparisons are secondary to the frozen full-profile question."}
    engine.dump(OUT/"score.json", result)
    lines = ["# New N400 homothetic three-shape stream", "",
        "Eight million new shared-counter replicas per shape pair, 400 aligned batches. Same three moduli and O, with all period entries doubled relative to N100.", "",
        "| shape | A_top | E_top | C | W |", "|---|---:|---:|---:|---:|"]
    for i, shape in enumerate(contract["shapes"]):
        values = [f"{mean[4*i+j]:.9g} +/- {np.sqrt(covariance[4*i+j,4*i+j]):.3g}" for j in range(4)]
        lines.append("| "+shape["name"]+" | "+" | ".join(values)+" |")
    lines += ["", "## Fixed-p secondary shape comparisons", "",
              "| model | chi-square / df | nominal p |", "|---|---:|---:|"]
    for name, row in result["same_area_models"].items():
        lines.append(f"| {name} | {row['chi_square']:.8g} / {row['df']} | {row['p_value']:.6g} |")
    lines += ["", result["inference"], "", "The complete threshold histograms are the main scientific product: full-profile invariants and scale changes are scored separately. No new field identification follows from these secondary comparisons.", ""]
    (OUT/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
