#!/usr/bin/env python3
"""Put threshold clocks, rank probabilities and birth insertions in one basis.

The committed production archives contain two complementary sufficient views
of the same Newman--Ziff batches:

* marginal K1/K2 histograms, which reconstruct F1,F2 and their p derivatives;
* joint moments through K1*K2, which reconstruct the C/W clock plane.

Joining on (N, orientation, batch) restores the complete batch-estimator
covariance between these views.  The event-level nonlinear cross-covariance
still requires the joint K1/K2 histogram and is reported as unavailable.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
from statistics import median
from typing import Any, Iterable, Sequence

import mpmath as mp

from analyze_matching_parity_derivatives_fast import (
    H,
    combine,
    cos4,
    read as read_histograms,
    tail_derivative,
    tail_expectation,
)


ROOT = Path(__file__).resolve().parents[1]

DEFAULT_ARCHIVES = (
    (
        "P49-N130",
        "results/server-20260828/P49-fullcurve-doubling-100m/raw/n130.hist.csv",
        "results/server-20260828/P49-fullcurve-doubling-100m/raw/n130.moments.csv",
    ),
    (
        "P50-N145",
        "results/server-20260829/P50-n145-n290-fullcurve/raw/n145_100m.hist.csv",
        "results/server-20260829/P50-n145-n290-fullcurve/raw/n145_100m.moments.csv",
    ),
    (
        "P49-N170",
        "results/server-20260828/P49-fullcurve-doubling-100m/raw/n170.hist.csv",
        "results/server-20260828/P49-fullcurve-doubling-100m/raw/n170.moments.csv",
    ),
    (
        "P43-N185",
        "results/server-20260828/P43-heldout-fullcurve-500m/raw/n185.hist.csv",
        "results/server-20260828/P43-heldout-fullcurve-500m/raw/n185.moments.csv",
    ),
    (
        "P154-N260",
        "results/server-20260829/P154-norm4-variance-pilot/raw/n260_10m.hist.csv",
        "results/server-20260829/P154-norm4-variance-pilot/raw/n260_10m.moments.csv",
    ),
    (
        "P43-N265",
        "results/server-20260828/P43-heldout-fullcurve-500m/raw/n265.hist.csv",
        "results/server-20260828/P43-heldout-fullcurve-500m/raw/n265.moments.csv",
    ),
    (
        "P50-N290",
        "results/server-20260829/P50-n145-n290-fullcurve/raw/n290_100m.hist.csv",
        "results/server-20260829/P50-n145-n290-fullcurve/raw/n290_100m.moments.csv",
    ),
    (
        "P57-N325",
        "results/server-20260829/P57-norm5-500m/raw/n325_500m.hist.csv",
        "results/server-20260829/P57-norm5-500m/raw/n325_500m.moments.csv",
    ),
    (
        "P154-N340",
        "results/server-20260829/P154-norm4-variance-pilot/raw/n340_10m.hist.csv",
        "results/server-20260829/P154-norm4-variance-pilot/raw/n340_10m.moments.csv",
    ),
    (
        "P57-N425",
        "results/server-20260829/P57-norm5-500m/raw/n425_500m.hist.csv",
        "results/server-20260829/P57-norm5-500m/raw/n425_500m.moments.csv",
    ),
)


METRICS = (
    "P4_C",
    "P4_W",
    "P4_A_top",
    "P4_E_top",
    "P4_S_birth",
    "P4_D_birth",
    "P4_F1",
    "P4_F2",
    "P4_f1",
    "P4_f2",
    "P4_S_historical",
    "P4_D_historical",
    "P4_S_historical_prime",
    "P4_D_historical_prime",
)


@dataclass(frozen=True)
class MomentRow:
    n: int
    a: int
    b: int
    orientation: str
    batch: int
    samples: int
    sum_k1: int
    sum_k2: int
    sum_k1_squared: int
    sum_k2_squared: int
    sum_product: int

    def clock(self) -> tuple[mp.mpf, mp.mpf]:
        scale = mp.mpf(self.n + 1) * self.samples
        c = mp.mpf(self.sum_k1 + self.sum_k2) / (2 * scale)
        w = mp.mpf(self.sum_k2 - self.sum_k1) / scale
        return c, w

    def clock_covariance(self) -> dict[str, Any]:
        c, w = self.clock()
        scale_squared = mp.mpf((self.n + 1) ** 2) * self.samples
        c_squared = mp.mpf(
            self.sum_k1_squared + 2 * self.sum_product + self.sum_k2_squared
        ) / (4 * scale_squared)
        w_squared = mp.mpf(
            self.sum_k1_squared - 2 * self.sum_product + self.sum_k2_squared
        ) / scale_squared
        c_times_w = mp.mpf(
            self.sum_k2_squared - self.sum_k1_squared
        ) / (2 * scale_squared)
        variance_c = c_squared - c * c
        variance_w = w_squared - w * w
        covariance = c_times_w - c * w
        return {
            "mean": {"C": _number(c), "W": _number(w)},
            "metric_order": ["C", "W"],
            "covariance": [
                [_number(variance_c), _number(covariance)],
                [_number(covariance), _number(variance_w)],
            ],
        }


def _number(value: mp.mpf | float | int, digits: int = 20) -> float:
    return float(mp.nstr(value, digits))


def read_moments(path: Path) -> dict[tuple[int, str, int], MomentRow]:
    rows: dict[tuple[int, str, int], MomentRow] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            row = MomentRow(
                n=int(raw["n"]),
                a=int(raw["a"]),
                b=int(raw["b"]),
                orientation=raw["orientation"],
                batch=int(raw["batch"]),
                samples=int(raw["samples"]),
                sum_k1=int(raw["sum_kminus"]),
                sum_k2=int(raw["sum_kplus"]),
                sum_k1_squared=int(raw["sum_kminus2"]),
                sum_k2_squared=int(raw["sum_kplus2"]),
                sum_product=int(raw["sum_product"]),
            )
            key = (row.n, row.orientation, row.batch)
            if key in rows:
                raise ValueError(f"duplicate moment batch {key} in {path}")
            rows[key] = row
    if not rows:
        raise ValueError(f"empty moments file: {path}")
    return rows


def combine_moments(rows: Sequence[MomentRow]) -> MomentRow:
    if not rows:
        raise ValueError("cannot combine an empty moments sequence")
    first = rows[0]
    if any(
        (row.n, row.a, row.b, row.orientation)
        != (first.n, first.a, first.b, first.orientation)
        for row in rows
    ):
        raise ValueError("incompatible moment rows")
    return MomentRow(
        n=first.n,
        a=first.a,
        b=first.b,
        orientation=first.orientation,
        batch=-1,
        samples=sum(row.samples for row in rows),
        sum_k1=sum(row.sum_k1 for row in rows),
        sum_k2=sum(row.sum_k2 for row in rows),
        sum_k1_squared=sum(row.sum_k1_squared for row in rows),
        sum_k2_squared=sum(row.sum_k2_squared for row in rows),
        sum_product=sum(row.sum_product for row in rows),
    )


def threshold_second_derivative(
    n: int, counts: Sequence[int], samples: int, p: mp.mpf
) -> mp.mpf:
    """Derivative of the beta-density mixture for one threshold endpoint."""

    if not 0 < p < 1:
        return mp.mpf(0)
    q = 1 - p
    density = n * q ** (n - 1)
    total = mp.mpf(0)
    for rank in range(1, n + 1):
        logarithmic_derivative = (rank - 1) / p - (n - rank) / q
        total += counts[rank] * density * logarithmic_derivative
        if rank < n:
            density *= mp.mpf(n - rank) * p / (rank * q)
    return total / samples


def endpoint_observables(row: H, p: mp.mpf) -> dict[str, mp.mpf]:
    f1_cdf = tail_expectation(row.n, row.minus, row.samples, p)
    f2_cdf = tail_expectation(row.n, row.plus, row.samples, p)
    f1_pdf = tail_derivative(row.n, row.minus, row.samples, p)
    f2_pdf = tail_derivative(row.n, row.plus, row.samples, p)
    f1_second = threshold_second_derivative(
        row.n, row.minus, row.samples, p
    )
    f2_second = threshold_second_derivative(
        row.n, row.plus, row.samples, p
    )
    a_top = f1_cdf + f2_cdf - 1
    e_top = 1 - f1_cdf + f2_cdf
    s_birth = f1_pdf + f2_pdf
    d_birth = f2_pdf - f1_pdf
    return {
        "F1": f1_cdf,
        "F2": f2_cdf,
        "f1": f1_pdf,
        "f2": f2_pdf,
        "A_top": a_top,
        "E_top": e_top,
        "S_birth": s_birth,
        "D_birth": d_birth,
        "S_birth_prime": f1_second + f2_second,
        "D_birth_prime": f2_second - f1_second,
        # Historical notation from analyze_matching_parity_derivatives.py.
        "S_historical": e_top / 2,
        "D_historical": a_top / 2,
        "S_historical_prime": d_birth / 2,
        "D_historical_prime": s_birth / 2,
    }


def intrinsic_center(first: H, second: H, iterations: int = 90) -> mp.mpf:
    lower = mp.mpf(0)
    upper = mp.mpf(1)
    for _ in range(iterations):
        p = (lower + upper) / 2
        average_a = (
            endpoint_observables(first, p)["A_top"]
            + endpoint_observables(second, p)["A_top"]
        ) / 2
        if average_a < 0:
            lower = p
        else:
            upper = p
    return (lower + upper) / 2


def _projected_state(
    first: dict[str, mp.mpf], second: dict[str, mp.mpf], delta_cos4: mp.mpf
) -> dict[str, mp.mpf]:
    def contrast(name: str) -> mp.mpf:
        return (first[name] - second[name]) / delta_cos4

    return {
        "P4_A_top": contrast("A_top"),
        "P4_E_top": contrast("E_top"),
        "P4_S_birth": contrast("S_birth"),
        "P4_D_birth": contrast("D_birth"),
        "P4_F1": contrast("F1"),
        "P4_F2": contrast("F2"),
        "P4_f1": contrast("f1"),
        "P4_f2": contrast("f2"),
        "P4_S_historical": contrast("S_historical"),
        "P4_D_historical": contrast("D_historical"),
        "P4_S_historical_prime": contrast("S_historical_prime"),
        "P4_D_historical_prime": contrast("D_historical_prime"),
    }


def _projected_p_derivative(
    first: dict[str, mp.mpf], second: dict[str, mp.mpf], delta_cos4: mp.mpf
) -> dict[str, mp.mpf]:
    def contrast(left: mp.mpf, right: mp.mpf) -> mp.mpf:
        return (left - right) / delta_cos4

    return {
        "P4_C": mp.mpf(0),
        "P4_W": mp.mpf(0),
        "P4_A_top": contrast(first["S_birth"], second["S_birth"]),
        "P4_E_top": contrast(first["D_birth"], second["D_birth"]),
        "P4_S_birth": contrast(
            first["S_birth_prime"], second["S_birth_prime"]
        ),
        "P4_D_birth": contrast(
            first["D_birth_prime"], second["D_birth_prime"]
        ),
        "P4_F1": contrast(first["f1"], second["f1"]),
        "P4_F2": contrast(first["f2"], second["f2"]),
        "P4_f1": contrast(
            (first["S_birth_prime"] - first["D_birth_prime"]) / 2,
            (second["S_birth_prime"] - second["D_birth_prime"]) / 2,
        ),
        "P4_f2": contrast(
            (first["S_birth_prime"] + first["D_birth_prime"]) / 2,
            (second["S_birth_prime"] + second["D_birth_prime"]) / 2,
        ),
        "P4_S_historical": contrast(
            first["D_birth"] / 2, second["D_birth"] / 2
        ),
        "P4_D_historical": contrast(
            first["S_birth"] / 2, second["S_birth"] / 2
        ),
        "P4_S_historical_prime": contrast(
            first["D_birth_prime"] / 2,
            second["D_birth_prime"] / 2,
        ),
        "P4_D_historical_prime": contrast(
            first["S_birth_prime"] / 2,
            second["S_birth_prime"] / 2,
        ),
    }


def _clock_projection(
    first: MomentRow, second: MomentRow, delta_cos4: mp.mpf
) -> dict[str, mp.mpf]:
    first_c, first_w = first.clock()
    second_c, second_w = second.clock()
    return {
        "P4_C": (first_c - second_c) / delta_cos4,
        "P4_W": (first_w - second_w) / delta_cos4,
    }


def _mean(values: Iterable[mp.mpf]) -> mp.mpf:
    rows = list(values)
    return mp.fsum(rows) / len(rows)


def _covariance_of_mean(influences: Sequence[Sequence[mp.mpf]]) -> list[list[mp.mpf]]:
    batches = len(influences)
    if batches < 2:
        raise ValueError("at least two common-field batches are required")
    return [
        [
            mp.fsum(row[i] * row[j] for row in influences)
            / (batches * (batches - 1))
            for j in range(len(METRICS))
        ]
        for i in range(len(METRICS))
    ]


def _matrix_payload(matrix: Sequence[Sequence[mp.mpf]]) -> list[list[float]]:
    return [[_number(value, 16) for value in row] for row in matrix]


def _identity_residuals(point: dict[str, mp.mpf]) -> dict[str, mp.mpf]:
    return {
        "A_equals_2D_historical": point["P4_A_top"]
        - 2 * point["P4_D_historical"],
        "E_equals_2S_historical": point["P4_E_top"]
        - 2 * point["P4_S_historical"],
        "S_birth_equals_2D_historical_prime": point["P4_S_birth"]
        - 2 * point["P4_D_historical_prime"],
        "D_birth_equals_2S_historical_prime": point["P4_D_birth"]
        - 2 * point["P4_S_historical_prime"],
        "F1_equals_half_A_minus_E": point["P4_F1"]
        - (point["P4_A_top"] - point["P4_E_top"]) / 2,
        "F2_equals_half_A_plus_E": point["P4_F2"]
        - (point["P4_A_top"] + point["P4_E_top"]) / 2,
        "f1_equals_half_S_minus_D": point["P4_f1"]
        - (point["P4_S_birth"] - point["P4_D_birth"]) / 2,
        "f2_equals_half_S_plus_D": point["P4_f2"]
        - (point["P4_S_birth"] + point["P4_D_birth"]) / 2,
    }


def analyze_archive(
    archive_id: str, histogram_path: Path, moment_path: Path
) -> dict[str, Any]:
    histograms = read_histograms(histogram_path)
    moments = read_moments(moment_path)
    if set(histograms) != set(moments):
        missing_hist = sorted(set(moments) - set(histograms))
        missing_moments = sorted(set(histograms) - set(moments))
        raise ValueError(
            f"unaligned archives for {archive_id}: "
            f"missing_hist={missing_hist[:3]}, missing_moments={missing_moments[:3]}"
        )
    sizes = {key[0] for key in histograms}
    if len(sizes) != 1:
        raise ValueError(f"{archive_id} contains multiple sizes")
    n = sizes.pop()
    grouped_hist = {
        orientation: [
            histograms[key]
            for key in sorted(histograms)
            if key[1] == orientation
        ]
        for orientation in ("first", "second")
    }
    batches = sorted(
        set(key[2] for key in histograms if key[1] == "first")
        & set(key[2] for key in histograms if key[1] == "second")
    )
    if not batches or any(len(grouped_hist[side]) != len(batches) for side in grouped_hist):
        raise ValueError(f"unaligned orientation batches in {archive_id}")
    first_total = combine(grouped_hist["first"])
    second_total = combine(grouped_hist["second"])
    p0 = intrinsic_center(first_total, second_total)
    delta_cos4 = cos4(first_total.a, first_total.b) - cos4(
        second_total.a, second_total.b
    )
    if delta_cos4 == 0:
        raise ValueError(f"zero H4 leverage in {archive_id}")

    total_first = endpoint_observables(first_total, p0)
    total_second = endpoint_observables(second_total, p0)
    point = _projected_state(total_first, total_second, delta_cos4)
    point.update(
        {
            name: _mean(
                _clock_projection(
                    moments[(n, "first", batch)],
                    moments[(n, "second", batch)],
                    delta_cos4,
                )[name]
                for batch in batches
            )
            for name in ("P4_C", "P4_W")
        }
    )
    p_derivative = _projected_p_derivative(
        total_first, total_second, delta_cos4
    )
    center_slope = (total_first["S_birth"] + total_second["S_birth"]) / 2
    if center_slope <= 0:
        raise ValueError("intrinsic matching center is not increasing")

    fixed_influences: list[list[mp.mpf]] = []
    center_influences: list[list[mp.mpf]] = []
    center_shifts: list[mp.mpf] = []
    for batch in batches:
        first_hist = histograms[(n, "first", batch)]
        second_hist = histograms[(n, "second", batch)]
        first_obs = endpoint_observables(first_hist, p0)
        second_obs = endpoint_observables(second_hist, p0)
        batch_point = _projected_state(first_obs, second_obs, delta_cos4)
        batch_point.update(
            _clock_projection(
                moments[(n, "first", batch)],
                moments[(n, "second", batch)],
                delta_cos4,
            )
        )
        fixed = [batch_point[name] - point[name] for name in METRICS]
        center_equation = (first_obs["A_top"] + second_obs["A_top"]) / 2
        center_shift = -center_equation / center_slope
        adjusted = [
            value + p_derivative[name] * center_shift
            for name, value in zip(METRICS, fixed)
        ]
        fixed_influences.append(fixed)
        center_influences.append(adjusted)
        center_shifts.append(center_shift)

    fixed_covariance = _covariance_of_mean(fixed_influences)
    center_covariance = _covariance_of_mean(center_influences)
    standard_errors = {
        name: mp.sqrt(max(mp.mpf(0), center_covariance[index][index]))
        for index, name in enumerate(METRICS)
    }
    identity_residuals = _identity_residuals(point)
    if any(abs(value) > mp.mpf("1e-25") for value in identity_residuals.values()):
        raise ArithmeticError(f"basis identity failed in {archive_id}")

    state_denominator = abs(point["P4_A_top"]) + abs(point["P4_E_top"])
    clock_denominator = abs(point["P4_C"]) + abs(point["P4_W"]) / 2
    density_denominator = abs(point["P4_D_birth"]) + abs(point["P4_S_birth"])
    total_moments = {
        orientation: combine_moments(
            [moments[(n, orientation, batch)] for batch in batches]
        )
        for orientation in ("first", "second")
    }
    return {
        "id": archive_id,
        "N": n,
        "histogram": str(histogram_path.relative_to(ROOT)),
        "moments": str(moment_path.relative_to(ROOT)),
        "orientations": {
            "first": [first_total.a, first_total.b],
            "second": [second_total.a, second_total.b],
        },
        "batches": len(batches),
        "samples_per_orientation": first_total.samples,
        "p0": _number(p0),
        "delta_cos4": _number(delta_cos4),
        "point": {name: _number(point[name]) for name in METRICS},
        "standard_error_center_influence": {
            name: _number(standard_errors[name], 16) for name in METRICS
        },
        "event_clock_covariance_by_orientation": {
            orientation: row.clock_covariance()
            for orientation, row in total_moments.items()
        },
        "covariance_metric_order": list(METRICS),
        "covariance_fixed_center_exact_batch_estimator": _matrix_payload(
            fixed_covariance
        ),
        "covariance_intrinsic_center_first_order_influence": _matrix_payload(
            center_covariance
        ),
        "center_influence": {
            "definition": "IF_p(batch)=-g_batch(p0)/mean_g_prime(p0)",
            "standard_error": _number(
                mp.sqrt(
                    mp.fsum(value * value for value in center_shifts)
                    / (len(center_shifts) * (len(center_shifts) - 1))
                ),
                16,
            ),
        },
        "exact_basis_identity_residuals": {
            name: _number(value) for name, value in identity_residuals.items()
        },
        "parity_direction_diagnostics": {
            "clock_odd_share_C_vs_W_over_2": _number(
                abs(point["P4_C"]) / clock_denominator
                if clock_denominator
                else mp.nan
            ),
            "state_odd_share_A_vs_E": _number(
                abs(point["P4_A_top"]) / state_denominator
                if state_denominator
                else mp.nan
            ),
            "density_odd_share_D_vs_S": _number(
                abs(point["P4_D_birth"]) / density_denominator
                if density_denominator
                else mp.nan
            ),
            "K2_state_cancellation_fraction": _number(
                abs(point["P4_F2"]) / (state_denominator / 2)
                if state_denominator
                else mp.nan
            ),
            "K2_clock_cancellation_fraction": _number(
                abs(point["P4_C"] + point["P4_W"] / 2) / clock_denominator
                if clock_denominator
                else mp.nan
            ),
            "K2_density_cancellation_fraction": _number(
                abs(point["P4_f2"]) / (density_denominator / 2)
                if density_denominator
                else mp.nan
            ),
        },
    }


def _sign_count(records: Sequence[dict[str, Any]], metric: str) -> dict[str, int]:
    values = [record["point"][metric] for record in records]
    return {
        "negative": sum(value < 0 for value in values),
        "zero": sum(value == 0 for value in values),
        "positive": sum(value > 0 for value in values),
        "total": len(values),
    }


def _maximin_geometry(records: Sequence[dict[str, Any]]) -> dict[str, Any]:
    rows = []
    for record in records:
        leverage_per_sqrt_site_cost = abs(record["delta_cos4"]) / math.sqrt(
            record["N"]
        )
        rows.append(
            {
                "id": record["id"],
                "N": record["N"],
                "orientations": record["orientations"],
                "delta_cos4": record["delta_cos4"],
                "deterministic_H4_leverage_per_sqrt_site_cost": leverage_per_sqrt_site_cost,
            }
        )
    rows.sort(
        key=lambda row: row["deterministic_H4_leverage_per_sqrt_site_cost"],
        reverse=True,
    )
    return {
        "single_size_mark_acquisition": {
            "selected": rows[0],
            "criterion": "maximize |Delta cos4|/sqrt(N) over committed production geometries, without using a new marked outcome",
            "ranking": rows,
        },
        "minimal_radial_campaign": {
            "selected": "N65_to_N130_q2",
            "parent": {"N": 65, "first": [8, 1], "second": [7, 4]},
            "child": {"N": 130, "first": [9, 7], "second": [11, 3]},
            "multiplier": [1, -1],
            "exact_line_character_transfer": "chi4((1-i) ell)=-chi4(ell)",
            "reason": "lowest-cost frozen radial edge with exact H4 sign reversal; acquire the marked tuple on both sizes to distinguish translation, lifetime and line polarization",
        },
    }


def build_report() -> dict[str, Any]:
    # Keep checked artifacts deterministic when build_report is imported by tests.
    mp.mp.dps = max(mp.mp.dps, 35)
    records = [
        analyze_archive(archive_id, ROOT / hist, ROOT / moments)
        for archive_id, hist, moments in DEFAULT_ARCHIVES
    ]
    records.sort(key=lambda row: row["N"])
    high_statistics = [row for row in records if not row["id"].startswith("P154")]
    diagnostics = (
        "clock_odd_share_C_vs_W_over_2",
        "state_odd_share_A_vs_E",
        "density_odd_share_D_vs_S",
        "K2_state_cancellation_fraction",
        "K2_clock_cancellation_fraction",
        "K2_density_cancellation_fraction",
    )
    medians = {
        name: median(
            row["parity_direction_diagnostics"][name] for row in high_statistics
        )
        for name in diagnostics
    }
    return {
        "schema": "matching-one/rank-plane-crosswalk/v1",
        "issues": [28, 43, 48, 215, 269, 275, 276],
        "status": "exact_basis_plus_same_stream_archive_reconstruction",
        "exact_crosswalk": {
            "linear_maps_for_orientation_contrasts": {
                "F1_F2_to_P0_P1_P2": {
                    "row_order": ["P0", "P1", "P2"],
                    "column_order": ["F1", "F2"],
                    "matrix": [[-1, 0], [1, -1], [0, 1]],
                },
                "F1_F2_to_A_E": {
                    "row_order": ["A_top", "E_top"],
                    "column_order": ["F1", "F2"],
                    "matrix": [[1, 1], [-1, 1]],
                    "affine_offsets_before_contrast": [-1, 1],
                },
                "f01_f12_to_D_S": {
                    "row_order": ["D_birth", "S_birth"],
                    "column_order": ["f01", "f12"],
                    "matrix": [[-1, 1], [1, 1]],
                },
                "K1_K2_to_c_W": {
                    "row_order": ["c", "W"],
                    "column_order": ["K1", "K2"],
                    "matrix_numerator": [[0.5, 0.5], [-1, 1]],
                    "common_denominator": "N+1",
                    "affine_offsets_before_contrast": [-0.5, 0],
                },
                "canonical_to_historical": {
                    "relations": [
                        "S_historical=E_top/2",
                        "D_historical=A_top/2",
                        "S_historical_prime=D_birth/2",
                        "D_historical_prime=S_birth/2",
                    ]
                },
            },
            "threshold_state": {
                "definitions": "F1=P(K1<=Bin(N,p)); F2=P(K2<=Bin(N,p)); P0=1-F1; P1=F1-F2; P2=F2",
                "canonical": "A_top=P2-P0=F1+F2-1; E_top=P0+P2=1-F1+F2",
                "historical_notation": "S_historical=E_top/2; D_historical=A_top/2",
            },
            "birth_density": {
                "definitions": "f01=F1_prime; f12=F2_prime",
                "canonical": "S_birth=f01+f12=A_top_prime; D_birth=f12-f01=E_top_prime=-P1_prime",
                "historical_notation": "D_historical_prime=S_birth/2; S_historical_prime=D_birth/2",
            },
            "clock": {
                "definitions": "c=C-1/2=(K1+K2-(N+1))/(2(N+1)); W=(K2-K1)/(N+1)",
                "inverse": "K1/(N+1)=C-W/2; K2/(N+1)=C+W/2",
            },
            "complement_parity": {
                "odd": ["c", "A_top", "D_birth"],
                "even": ["W", "E_top", "S_birth"],
                "endpoint_mixing": "F2=(A_top+E_top)/2 and f12=(S_birth+D_birth)/2 are not parity eigenvectors",
            },
            "marked_extension": {
                "strict_plateau": "for K1<K2, one-hot sectors (ell,iota) split W=sum W_(ell,iota) and D_birth=sum D_(ell,iota)",
                "spin4": "chi4(ell)W is complement-even; chi4(ell)D_birth is complement-odd",
                "direct_0_to_2": "K1=K2 has W=0 and D_birth=0 but contributes to an additional unpolarized S_birth_direct sector",
                "iota_boundary": "iota preservation is currently a tiny-oracle fact, not a general theorem",
            },
        },
        "archive_recoverability": {
            "exact": [
                "C/W means and event-level 2x2 covariance from K1,K2 first and second joint moments",
                "F1/F2 and all p derivatives from marginal endpoint histograms",
                "full fixed-center covariance of batch estimators across C,W,A,E,S_birth,D_birth by joining common batch ids",
            ],
            "first_order": [
                "full intrinsic-center covariance after the displayed M-estimator center influence correction"
            ],
            "not_recoverable": [
                "event-level nonlinear covariance between C/W and endpoint CDF or density kernels in production archives without joint K1/K2 histograms",
                "any ell/iota-resolved mean or covariance",
                "chi4(ell) S/D and same-sample A_top times J_D4 connected moments",
            ],
        },
        "datasets": records,
        "high_statistics_summary": {
            "selection": "exclude only P154 10M variance pilots",
            "sizes": [row["N"] for row in high_statistics],
            "sign_counts": {
                metric: _sign_count(high_statistics, metric)
                for metric in (
                    "P4_C",
                    "P4_W",
                    "P4_A_top",
                    "P4_E_top",
                    "P4_S_birth",
                    "P4_D_birth",
                    "P4_F2",
                )
            },
            "median_direction_diagnostics": medians,
            "K2_cancellation_answer": "No. K2/F2 is the mixed-parity vector (A_top+E_top)/2. Its historical cancellation is a cross-parity state-level cancellation, not the complement-even W/S_birth projector.",
            "H4_direction_answer": "The clock and density layers are stably odd-leading in C and D_birth, but W and S_birth are nonzero companion directions. The rank plane is two-component, not a pure Alexander-odd line.",
        },
        "maximin_next_acquisition": _maximin_geometry(high_statistics),
        "minimal_sufficient_statistics": {
            "threshold_stream": {
                "one_sparse_row": [
                    "lineage",
                    "orientation",
                    "batch",
                    "K1",
                    "K2",
                    "ell_u",
                    "ell_v",
                    "iota",
                    "count",
                ],
                "null_rule": "ell/iota are null exactly for direct K1=K2 activations",
                "reconstructs": "event-level covariance of C,W and every p-resolved unmarked/chi4-marked S_birth,D_birth kernel",
            },
            "connected_source_extension": [
                "sum q",
                "sum q^2",
                "sum J_D4_re/im",
                "sum q*J_D4_re/im",
                "sum J_S4_re/im as the complement sign control",
            ],
        },
        "scientific_card": [
            "MECHANISM SPACE: K2 cancellation is resolved as A_top/E_top cross-parity interference; the exact eigenplanes are (C,W), (A,E), and (D_birth,S_birth).",
            "NOT PROVED: archives do not establish a one-dimensional Q4 field, an ell/iota continuum overlap, or an asymptotic exponent for W/S_birth.",
            "OBSERVER-SECTOR-SOURCE-GEOMETRY: A_top/E_top | odd/even Alexander sectors | chi4(ell)D/S sources | Gaussian same-N orientation pairs.",
            "DEPENDENCY GROUP: all reconstructed contrasts and covariances are views of the same threshold-rank streams, not independent evidence rows.",
            "UPWEIGHT OBSERVATION: acquire sparse joint (K1,K2,ell,iota) on the q2 N65->N130 sign-flip edge, with qJ_D4 cross moments when the connected #275 coupling is targeted.",
        ],
    }


def render_markdown(report: dict[str, Any]) -> str:
    summary = report["high_statistics_summary"]
    medians = summary["median_direction_diagnostics"]
    selected = report["maximin_next_acquisition"]["single_size_mark_acquisition"][
        "selected"
    ]
    lines = [
        "# Rank-plane crosswalk: clocks, projectors and typed insertions",
        "",
        "## Exact basis",
        "",
        "```text",
        "clock odd/even:       c=C-1/2                 W=(K2-K1)/(N+1)",
        "state odd/even:       A=P2-P0                 E=P0+P2=1-P1",
        "birth odd/even:       D=f12-f01=E'            S=f01+f12=A'",
        "historical notation:  D_hist=A/2, S_hist=E/2, S_hist'=D/2, D_hist'=S/2",
        "endpoint mixtures:    F2=(A+E)/2              f12=(S+D)/2",
        "```",
        "",
        "Complement/reversal makes `c,A,D` odd and `W,E,S` even. Therefore",
        "the historical `K2`/`F2` cancellation cannot be the even `W/S`",
        "projector: `F2` is exactly a mixed-parity endpoint.",
        "",
        "## Archive result",
        "",
        f"The joined high-statistics streams cover `N={summary['sizes']}`.",
        f"The median odd shares are `{medians['clock_odd_share_C_vs_W_over_2']:.3f}`",
        f"in the clock plane, `{medians['state_odd_share_A_vs_E']:.3f}` in the",
        f"state plane, and `{medians['density_odd_share_D_vs_S']:.3f}` in the",
        "birth-density plane. The median K2 cancellation fraction is",
        f"`{medians['K2_state_cancellation_fraction']:.3f}` at the state layer,",
        f"versus `{medians['K2_clock_cancellation_fraction']:.3f}` for its clock",
        f"endpoint and `{medians['K2_density_cancellation_fraction']:.3f}` for its",
        "density endpoint. The strong cancellation is therefore localized to",
        "the state value at the intrinsic center; it is not a universal K2 null.",
        "",
        "Each dataset in the JSON carries a complete `14 x 14` same-batch",
        "covariance in this basis. The fixed-center matrix is an exact covariance",
        "of the archived batch estimators. The intrinsic-center matrix adds the",
        "displayed first-order M-estimator influence of the fitted center.",
        "",
        "## What is and is not reconstructible",
        "",
        "The production moments recover event-level `Cov(C,W)`, while marginal",
        "histograms recover all endpoint CDF/density means. Common batch IDs",
        "recover their estimator covariance. Event-level nonlinear cross moments",
        "need the joint `K1,K2` histogram; `ell`, `iota`, `chi4(ell)D`, and `qJ`",
        "cannot be recovered from the old production files.",
        "",
        "## Next acquisition",
        "",
        f"Among committed single-size geometries, `{selected['id']}`",
        f"(`N={selected['N']}`, {selected['orientations']}) maximizes the",
        "outcome-free `|Delta cos4|/sqrt(N)` mark-acquisition proxy. The more",
        "informative minimal radial campaign is the frozen q2 edge",
        "`N65 -> N130`: its `(1-i)` multiplier flips `chi4(ell)` exactly and is",
        "the cheapest way to separate clock translation, rank-one lifetime and",
        "line polarization across scale.",
        "",
        "Store a sparse joint histogram of",
        "`(K1,K2,ell_u,ell_v,iota,count)` per aligned batch/orientation. Add",
        "`q`, `J_D4`, and `q*J_D4` batch sums only when the connected #275",
        "coupling is desired.",
        "",
        "## Scientific card",
        "",
    ]
    lines.extend(f"{index}. {line}" for index, line in enumerate(report["scientific_card"], 1))
    lines.extend(
        [
            "",
            "## Claim boundary",
            "",
            "The basis and parity statements are exact. Archive covariance is",
            "same-stream reuse; the intrinsic-center extension is first-order.",
            "Sign coherence and odd-leading language are mechanism diagnostics,",
            "not independent evidence or an asymptotic field identification.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--dps", type=int, default=35)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    report = build_report()
    rendered = (
        json.dumps(report, indent=2, sort_keys=True) + "\n"
        if args.format == "json"
        else render_markdown(report)
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
