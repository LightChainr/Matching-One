#!/usr/bin/env python3
"""Build the equal-area rectangle design for the P321 aspect-ratio bridge."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path


BASE_RECTANGLES = (
    ("rho_1", 12, 12, "primary"),
    ("rho_16_9", 9, 16, "primary"),
    ("rho_9_4", 8, 18, "primary"),
    ("rho_4", 6, 24, "primary"),
    ("rho_9", 4, 36, "endpoint_diagnostic"),
)
SCALES = (1, 2, 3)


def build_design() -> dict[str, object]:
    rows = []
    for scale in SCALES:
        expected_area = 144 * scale * scale
        for identifier, width0, height0, role in BASE_RECTANGLES:
            width, height = scale * width0, scale * height0
            if width * height != expected_area:
                raise AssertionError("rectangle family lost equal area")
            aspect = Fraction(height, width)
            rows.append(
                {
                    "id": f"s{scale}_{identifier}",
                    "scale": scale,
                    "N": expected_area,
                    "width": width,
                    "height": height,
                    "aspect_ratio": f"{aspect.numerator}/{aspect.denominator}",
                    "aspect_ratio_decimal": float(aspect),
                    "period_matrix_row_major": [[width, 0], [0, height]],
                    "role": role,
                }
            )
    return {
        "schema": "matching-one.p321-equal-area-rectangle-design.v1",
        "status": "geometry_frozen; variance pilot and sample count not authorized",
        "source_theory_commit": "81818ec",
        "observable": "matching root of P2-P0",
        "scaling_contract": {
            "leading_root_shift": "N^-2 at fixed aspect ratio (equivalently linear-scale^-4)",
            "first_relative_correction": "N^-1, producing root term N^-3",
            "width_amplitude_conversion": (
                "If p-pc=C_N(rho)*N^-2=C_width(rho)*n_width^-4 and "
                "N=rho*n_width^2, then C_width(rho)=C_N(rho)/rho^2. "
                "Apply this conversion before the rho->infinity TL comparison."
            ),
            "free_exponent_fit": False,
        },
        "common_randomness": {
            "rule": (
                "At each N, run the square as first matrix against every non-square "
                "second matrix with identical seed, counter interval, and batch boundaries."
            ),
            "exact_gate": (
                "The repeated square threshold histograms must be byte-identical; align "
                "delete-one-batch root vectors across all aspect ratios for full covariance."
            ),
        },
        "rows": rows,
        "primary_aspect_ratios": ["1", "16/9", "9/4", "4"],
        "endpoint_diagnostic": "rho=9 is not treated as the TL endpoint; it only constrains approach",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    design = build_design()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(design, indent=2) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
