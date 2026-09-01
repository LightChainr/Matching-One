#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


HERE = Path(__file__).resolve().parent
N = 25
DELTA = F(1152, 625)
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


@dataclass(frozen=True)
class Interval:
    lo: F
    hi: F

    @classmethod
    def of(cls, value):
        return value if isinstance(value, cls) else cls(F(value), F(value))

    def __add__(self, other):
        other = self.of(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -self.of(other)

    def __rsub__(self, other):
        return self.of(other) - self

    def __mul__(self, other):
        other = self.of(other)
        products = (
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        )
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.of(other)
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError("interval denominator contains zero")
        return self * Interval(1 / other.hi, 1 / other.lo)

    def __pow__(self, exponent):
        if exponent < 0:
            raise ValueError("nonnegative exponents only")
        value = Interval.of(1)
        for _ in range(exponent):
            value *= self
        return value


def interval_record(value: Interval):
    # Compact exact outward enclosure.  The underlying arithmetic remains
    # Fraction interval arithmetic; only serialization is rounded outward.
    scale = 10 ** 50
    lower = F(math.floor(value.lo * scale), scale)
    upper = F(math.ceil(value.hi * scale), scale)
    midpoint = (value.lo + value.hi) / 2
    return {
        "lower": str(lower),
        "upper": str(upper),
        "midpoint": float(midpoint),
        "width": float(upper - lower),
        "excludes_zero": lower > 0 or upper < 0,
        "serialization": "exact_outward_rounding_to_1e_minus50_rational_grid",
    }


def read_csv(path):
    with path.open(newline="") as handle:
        return [{key: int(value) for key, value in row.items()} for row in csv.DictReader(handle)]


def baseline_coefficients(path):
    result = {name: [F(0)] * (N + 1) for name in ("z", "q", "e")}
    for row in read_csv(path):
        k, q, count = row["k"], row["q"], row["count"]
        result["z"][k] += count
        result["q"][k] += count * q
        result["e"][k] += count * q * q
    assert result["z"] == [math.comb(N, k) for k in range(N + 1)]
    return result


def source_coefficients(path):
    rows = read_csv(path)
    assert len(rows) == N + 1
    result = {name: [F(0)] * (N + 1) for name in ("s", "qs", "es")}
    for k, row in enumerate(rows):
        assert row["k"] == k
        assert row["count"] == (math.comb(N - 1, k) if k < N else 0)
        for name, field in (("s", "sum_G16"), ("qs", "sum_G16_q"), ("es", "sum_G16_E")):
            result[name][k] = F(row[field], 16 * N)
    return result


def weight_derivatives(k, p):
    m = N - k
    pows = [p ** j for j in range(N + 1)]
    qs = [(1 - p) ** j for j in range(N + 1)]
    value = pows[k] * qs[m]
    first = (k * pows[k - 1] * qs[m] if k else 0) - (
        m * pows[k] * qs[m - 1] if m else 0
    )
    second = (
        (k * (k - 1) * pows[k - 2] * qs[m] if k >= 2 else 0)
        - (2 * k * m * pows[k - 1] * qs[m - 1] if k and m else 0)
        + (m * (m - 1) * pows[k] * qs[m - 2] if m >= 2 else 0)
    )
    return value, first, second


def jets(coefficients, p):
    output = {name: [Interval.of(0), Interval.of(0), Interval.of(0)] for name in coefficients}
    for k in range(N + 1):
        weights = weight_derivatives(k, p)
        for name, values in coefficients.items():
            for order in range(3):
                output[name][order] += values[k] * weights[order]
    return output


def geometry_packet(baseline_path, source_path, p):
    coefficients = baseline_coefficients(baseline_path)
    coefficients.update(source_coefficients(source_path))
    row = jets(coefficients, p)
    q, qp, qpp = row["q"]
    e, ep, epp = row["e"]
    s, sp, _ = row["s"]
    qs, qsp, _ = row["qs"]
    es, esp, _ = row["es"]
    return {
        "q": q,
        "q_p": qp,
        "q_pp": qpp,
        "e": e,
        "e_p": ep,
        "e_pp": epp,
        "s": s,
        "s_p": sp,
        "j_q": qs - q * s,
        "j_q_p": qsp - qp * s - q * sp,
        "j_e": es - e * s,
        "j_e_p": esp - ep * s - e * sp,
    }


def determinant(first, second):
    return first[0] * second[1] - first[1] * second[0]


def score(mode, p):
    packets = [
        geometry_packet(HERE / f"baseline-{name}.csv", HERE / f"{name}-{mode}.csv", p)
        for name in ("axis", "tilted")
    ]
    axis, tilted = packets
    thermal_rows = [axis["q_p"], axis["e_p"], tilted["q_p"], tilted["e_p"]]
    source_rows = [axis["j_q"], axis["j_e"], tilted["j_q"], tilted["j_e"]]
    labels = ["axis_q", "axis_E", "tilted_q", "tilted_E"]
    minors = []
    for i, j in combinations(range(4), 2):
        value = thermal_rows[i] * source_rows[j] - thermal_rows[j] * source_rows[i]
        minors.append({"rows": [labels[i], labels[j]], "determinant": interval_record(value)})

    M_p = (axis["q_p"] + tilted["q_p"]) / 2
    M_pp = (axis["q_pp"] + tilted["q_pp"]) / 2
    Y_p = (axis["e_p"] - tilted["e_p"]) / DELTA
    Y_pp = (axis["e_pp"] - tilted["e_pp"]) / DELTA
    jM = (axis["j_q"] + tilted["j_q"]) / 2
    jM_p = (axis["j_q_p"] + tilted["j_q_p"]) / 2
    jY = (axis["j_e"] - tilted["j_e"]) / DELTA
    jY_p = (axis["j_e_p"] - tilted["j_e_p"]) / DELTA
    projected_minor = M_p * jY - Y_p * jM
    fixed_root_source = projected_minor / M_p
    R = Y_p / M_p
    R_p = (Y_pp * M_p - Y_p * M_pp) / (M_p * M_p)
    T_p = jY_p - R * jM_p - R_p * jM
    mixed_hessian = T_p / M_p
    return {
        "mode": mode,
        "row_order": labels,
        "matrix": {
            "thermal_column": [interval_record(value) for value in thermal_rows],
            "source_column": [interval_record(value) for value in source_rows],
        },
        "all_six_row_minors": minors,
        "P4_root_projection": {
            "rows": ["M=(q_axis+q_tilted)/2", "Y=(E_axis-E_tilted)/Delta"],
            "matrix": [
                [interval_record(M_p), interval_record(jM)],
                [interval_record(Y_p), interval_record(jY)],
            ],
            "determinant": interval_record(projected_minor),
            "source_after_fixed_M_Schur_elimination": interval_record(fixed_root_source),
            "Schur_invariance": "source_column -> source_column-beta*thermal_column leaves every 2x2 minor unchanged",
        },
        "root_conditioned_mixed_hessian_Tp_over_Mp": interval_record(mixed_hessian),
        "nonzero_minor_count": sum(item["determinant"]["excludes_zero"] for item in minors),
        "all_row_minors_exclude_zero": all(item["determinant"]["excludes_zero"] for item in minors),
    }


def file_sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    root = json.loads((HERE / "baseline-root.json").read_text())
    p = Interval(
        F(root["root_enclosure"]["lower_fraction"]),
        F(root["root_enclosure"]["upper_fraction"]),
    )
    # Exact additive control on every K slice and every moment column.
    for geometry in ("axis", "tilted"):
        same = read_csv(HERE / f"{geometry}-clean_same.csv")
        reversed_rows = read_csv(HERE / f"{geometry}-clean_reversed.csv")
        total = read_csv(HERE / f"{geometry}-clean_total.csv")
        for left, right, combined in zip(same, reversed_rows, total):
            for field in ("sum_G16", "sum_G16_q", "sum_G16_E"):
                assert left[field] + right[field] == combined[field]
    modes = [score(mode, p) for mode in ("clean_same", "clean_reversed", "clean_total")]
    all_control = score("all", p)
    result = {
        "schema": "matching-one/p537-provisional-clean-two-bridge-minors/v1",
        "status": "exact_finite_counterexample_under_explicit_provisional_landing_contract",
        "landing_contract": {
            "name": "clean_two_bridge_six_block",
            "definition": "Bell-8 partition has exactly six blocks; exactly two blocks occur once at each marked four-port group; the other four blocks are singleton ports; the two shared landing positions are adjacent at both marks",
            "C4_orbits": {
                "clean_same": "cyclic order of the two shared labels agrees at the two marks",
                "clean_reversed": "cyclic order is reversed",
                "clean_total": "sum of both nonzero C4 orbits",
            },
            "kernel_check": "every retained physical record has g16=4; program aborts otherwise",
            "warning": "This is a minimal auditable interpretation introduced by this audit. The repository has not formally equated it with the phrase ordinary four-arm/no-extra-branch.",
        },
        "N": 25,
        "geometries": [[5, 0], [4, 3]],
        "root_p": interval_record(p),
        "Delta": str(DELTA),
        "modes": modes,
        "controls": {
            "configurations_per_geometry_per_mode": 2 ** 24,
            "six_exact_traversals": True,
            "clean_total_equals_orbit_sum_per_K_and_moment": True,
            "all_mode_axis_bitwise_matches_saved": (
                (HERE / "axis-all.csv").read_bytes() == (HERE / "original-axis.csv").read_bytes()
            ),
            "all_mode_tilted_bitwise_matches_saved": (
                (HERE / "tilted-all.csv").read_bytes() == (HERE / "original-tilted.csv").read_bytes()
            ),
            "all_mode_reproduced_full_J2_over_A": all_control[
                "root_conditioned_mixed_hessian_Tp_over_Mp"
            ],
            "published_full_J2_over_A_midpoint": -5.905706006949678e-05,
            "random_samples": 0,
            "cloud_jobs": 0,
            "artifact_is_self_contained": True,
        },
        "source_files": {
            path.name: file_sha256(path)
            for path in sorted(HERE.glob("*.csv"))
        },
    }
    (HERE / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    for row in modes:
        projection = row["P4_root_projection"]
        print(
            row["mode"],
            "nonzero minors", row["nonzero_minor_count"], "/6",
            "projected det", projection["determinant"]["midpoint"],
            "mixed", row["root_conditioned_mixed_hessian_Tp_over_Mp"]["midpoint"],
        )


if __name__ == "__main__":
    main()
