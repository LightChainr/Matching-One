#!/usr/bin/env python3
"""Score the exact P537 site-flip landing matrix from sufficient aggregates.

The input is one CSV with these columns:

    geometry,tau,alpha,source_component,k_minus,count,
    sum_q0,sum_q1,sum_e0,sum_e1,sum_a0,sum_a1,
    sum_q0a0,sum_q1a1,sum_e0a0,sum_e1a1

``tau=__GLOBAL__`` is reserved for the complete (not landing-filtered)
source profile for each ``geometry,source_component,k_minus``.  Those rows
determine ``mu_a[geometry,lambda]``, pooled ``jM[lambda]`` and
``beta[lambda]``.  All other tau values are disjoint landing cells; they may
refine a fixed component into fibre-dependent ``alpha`` labels, but Eq. (10)
is evaluated with the parent component's beta before cells are aggregated.
``notes/p537-finite-landing-transfer-definition.md``.

The a-fields are raw fixed-source-anchor alpha sums.  Translation fixes the
source endpoint 0, so every component is ``a^y=N^-1 g_0y``; this scorer
applies that one required 1/N factor to every field containing a.  It does
*not* divide again for the thermal site z.  By default the producer fixes one
NN z direction, and the scorer multiplies that representative by its four
element physical C4 orbit.  Pass ``--z-orbit-multiplicity 1`` only if all four
directions were already accumulated.  A source-absent cell must still be
written with zero a-fields and its nonzero count/q/E fields: its
``-beta_y*b_z`` Schur contribution is part of Eq. (10).

The canonical interface expects actual g values.  A producer retaining the
integer ``g16=16*g`` may pass those integer sums unchanged together with
``--a-raw-denominator 16``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

GLOBAL_TAU = "__GLOBAL__"
GEOMETRIES = ("axis", "tilted")
REQUIRED_FIELDS = {
    "geometry",
    "tau",
    "alpha",
    "source_component",
    "k_minus",
    "count",
    "sum_q0",
    "sum_q1",
    "sum_e0",
    "sum_e1",
    "sum_a0",
    "sum_a1",
    "sum_q0a0",
    "sum_q1a1",
    "sum_e0a0",
    "sum_e1a1",
}
A_FIELDS = {
    "sum_a0",
    "sum_a1",
    "sum_q0a0",
    "sum_q1a1",
    "sum_e0a0",
    "sum_e1a1",
}


@dataclass(frozen=True)
class Interval:
    lo: F
    hi: F

    @classmethod
    def of(cls, value: object) -> "Interval":
        return value if isinstance(value, cls) else cls(F(value), F(value))

    def __add__(self, other: object) -> "Interval":
        rhs = self.of(other)
        return Interval(self.lo + rhs.lo, self.hi + rhs.hi)

    __radd__ = __add__

    def __neg__(self) -> "Interval":
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other: object) -> "Interval":
        return self + -self.of(other)

    def __rsub__(self, other: object) -> "Interval":
        return self.of(other) - self

    def __mul__(self, other: object) -> "Interval":
        rhs = self.of(other)
        products = (
            self.lo * rhs.lo,
            self.lo * rhs.hi,
            self.hi * rhs.lo,
            self.hi * rhs.hi,
        )
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def __truediv__(self, other: object) -> "Interval":
        rhs = self.of(other)
        if rhs.lo <= 0 <= rhs.hi:
            raise ZeroDivisionError("interval denominator contains zero")
        return self * Interval(1 / rhs.hi, 1 / rhs.lo)

    def __pow__(self, exponent: int) -> "Interval":
        if exponent < 0:
            raise ValueError("nonnegative exponents only")
        result = Interval.of(1)
        for _ in range(exponent):
            result *= self
        return result


def parse_fraction(text: str) -> F:
    return F(text.strip())


def interval_record(value: Interval) -> dict[str, object]:
    scale = 10**50
    lower = F(math.floor(value.lo * scale), scale)
    upper = F(math.ceil(value.hi * scale), scale)
    midpoint = (value.lo + value.hi) / 2
    return {
        "lower": str(lower),
        "upper": str(upper),
        "midpoint": float(midpoint),
        "width": float(upper - lower),
        "excludes_zero": lower > 0 or upper < 0,
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_root(path: Path) -> Interval:
    payload = json.loads(path.read_text())
    root = payload.get("root_enclosure", payload.get("root_p"))
    if root is None:
        raise ValueError(f"{path}: missing root_enclosure/root_p")
    lower = root.get("lower_fraction", root.get("lower"))
    upper = root.get("upper_fraction", root.get("upper"))
    if lower is None or upper is None:
        raise ValueError(f"{path}: root enclosure lacks rational endpoints")
    return Interval(parse_fraction(lower), parse_fraction(upper))


def read_baseline(path: Path, n: int) -> dict[str, list[F]]:
    coefficients = {name: [F(0)] * (n + 1) for name in ("z", "q", "e")}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            k = int(row["k"])
            q = int(row["q"])
            count = int(row["count"])
            coefficients["z"][k] += count
            coefficients["q"][k] += count * q
            coefficients["e"][k] += count * q * q
    expected = [F(math.comb(n, k)) for k in range(n + 1)]
    if coefficients["z"] != expected:
        raise ValueError(f"{path}: baseline does not enumerate every N={n} state")
    return coefficients


def full_weight(k: int, n: int, p: Interval) -> Interval:
    return p**k * (1 - p) ** (n - k)


def offsite_weight(k: int, n: int, p: Interval) -> Interval:
    return p**k * (1 - p) ** (n - 1 - k)


def baseline_packet(coefficients: dict[str, list[F]], n: int, p: Interval) -> dict[str, Interval]:
    result = {name: Interval.of(0) for name in ("q", "e", "q_t", "e_t")}
    for k in range(n + 1):
        weight = full_weight(k, n, p)
        score = k - n * p
        result["q"] += coefficients["q"][k] * weight
        result["e"] += coefficients["e"][k] * weight
        result["q_t"] += coefficients["q"][k] * weight * score
        result["e_t"] += coefficients["e"][k] * weight * score
    return result


def read_aggregates(path: Path, n: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    seen: set[tuple[str, str, str, str, int]] = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or ())
        missing = REQUIRED_FIELDS - fields
        if missing:
            raise ValueError(f"{path}: missing fields {sorted(missing)}")
        for line_number, raw in enumerate(reader, start=2):
            geometry = raw["geometry"].strip()
            if geometry not in GEOMETRIES:
                raise ValueError(f"{path}:{line_number}: unknown geometry {geometry!r}")
            tau, alpha = raw["tau"].strip(), raw["alpha"].strip()
            component = raw["source_component"].strip()
            if not tau or not alpha or not component:
                raise ValueError(f"{path}:{line_number}: empty tau/alpha/source_component")
            k = int(raw["k_minus"])
            if not 0 <= k < n:
                raise ValueError(f"{path}:{line_number}: k_minus outside [0,{n - 1}]")
            key = (geometry, tau, alpha, component, k)
            if key in seen:
                raise ValueError(f"{path}:{line_number}: duplicate key {key}")
            seen.add(key)
            row: dict[str, object] = {
                "geometry": geometry,
                "tau": tau,
                "alpha": alpha,
                "component": component,
                "k": k,
            }
            for field in REQUIRED_FIELDS - {
                "geometry", "tau", "alpha", "source_component", "k_minus"
            }:
                row[field] = parse_fraction(raw[field])
            if row["count"] < 0:
                raise ValueError(f"{path}:{line_number}: negative count")
            rows.append(row)
    return rows


def allocated(row: dict[str, object], field: str, n: int, a_raw_denominator: int) -> F:
    value = row[field]
    assert isinstance(value, F)
    return value / (n * a_raw_denominator) if field in A_FIELDS else value


def eq10_state_term(
    weight_i: Interval,
    u_i: Interval,
    b_i: Interval,
    sum_h: Interval,
    sum_ha: Interval,
    beta_component: Interval,
    z_orbit_multiplicity: int,
) -> Interval:
    """One i=0/1 term, including the fixed-z physical orbit multiplier."""
    return z_orbit_multiplicity * weight_i * (
        u_i * sum_ha - beta_component * b_i * sum_h
    )


def source_global_packets(
    rows: list[dict[str, object]],
    baseline: dict[str, dict[str, Interval]],
    n: int,
    p: Interval,
    a_raw_denominator: int,
) -> tuple[dict[tuple[str, str], dict[str, Interval]], dict[str, Interval]]:
    global_rows = [row for row in rows if row["tau"] == GLOBAL_TAU]
    components = sorted({str(row["component"]) for row in global_rows})
    if not components:
        raise ValueError(f"no {GLOBAL_TAU!r} source rows")
    packets: dict[tuple[str, str], dict[str, Interval]] = {}
    for geometry in GEOMETRIES:
        for component in components:
            selected = [
                row for row in global_rows
                if row["geometry"] == geometry and row["component"] == component
            ]
            if not selected:
                raise ValueError(f"missing global source profile for {(geometry, component)}")
            mu_a = Interval.of(0)
            mean_qa = Interval.of(0)
            for row in selected:
                k = int(row["k"])
                weight = offsite_weight(k, n, p)
                a0 = allocated(row, "sum_a0", n, a_raw_denominator)
                a1 = allocated(row, "sum_a1", n, a_raw_denominator)
                qa0 = allocated(row, "sum_q0a0", n, a_raw_denominator)
                qa1 = allocated(row, "sum_q1a1", n, a_raw_denominator)
                mu_a += weight * ((1 - p) * a0 + p * a1)
                mean_qa += weight * ((1 - p) * qa0 + p * qa1)
            mu_q = baseline[geometry]["q"]
            packets[(geometry, component)] = {
                "mu_a": mu_a,
                "cov_q_a": mean_qa - mu_q * mu_a,
            }

    mt = (baseline["axis"]["q_t"] + baseline["tilted"]["q_t"]) / 2
    if mt.lo <= 0 <= mt.hi:
        raise ValueError("pooled matching slope interval contains zero")
    beta: dict[str, Interval] = {}
    for component in components:
        jm = (
            packets[("axis", component)]["cov_q_a"]
            + packets[("tilted", component)]["cov_q_a"]
        ) / 2
        beta[component] = jm / mt
        for geometry in GEOMETRIES:
            packets[(geometry, component)]["jM_component"] = jm
            packets[(geometry, component)]["beta_component"] = beta[component]
    return packets, beta


def score(args: argparse.Namespace) -> dict[str, object]:
    n = args.n
    if args.z_orbit_multiplicity <= 0:
        raise ValueError("z-orbit-multiplicity must be positive")
    if args.a_raw_denominator <= 0:
        raise ValueError("a-raw-denominator must be positive")
    delta = parse_fraction(args.delta)
    if delta == 0:
        raise ValueError("Delta must be nonzero")
    p = read_root(args.baseline_root)
    rows = read_aggregates(args.aggregates, n)
    coefficients = {
        "axis": read_baseline(args.baseline_axis, n),
        "tilted": read_baseline(args.baseline_tilted, n),
    }
    baseline = {
        geometry: baseline_packet(coefficients[geometry], n, p) for geometry in GEOMETRIES
    }
    mt = (baseline["axis"]["q_t"] + baseline["tilted"]["q_t"]) / 2
    yt = (baseline["axis"]["e_t"] - baseline["tilted"]["e_t"]) / delta
    r = yt / mt
    c = {"axis": F(1, 1) / delta, "tilted": -F(1, 1) / delta}
    mu_h = {
        geometry: 2 * c[geometry] * baseline[geometry]["e"]
        - r * baseline[geometry]["q"]
        for geometry in GEOMETRIES
    }
    source_packets, beta = source_global_packets(rows, baseline, n, p, args.a_raw_denominator)

    geometry_cells: dict[tuple[str, str, str], Interval] = defaultdict(lambda: Interval.of(0))
    for row in rows:
        tau = str(row["tau"])
        if tau == GLOBAL_TAU:
            continue
        geometry = str(row["geometry"])
        alpha = str(row["alpha"])
        component = str(row["component"])
        if (geometry, component) not in source_packets:
            raise ValueError(
                f"landing cell {(geometry, tau, alpha, component)} lacks global source profile"
            )
        k = int(row["k"])
        weight = offsite_weight(k, n, p)
        mu_a = source_packets[(geometry, component)]["mu_a"]
        beta_component = beta[component]
        s_minus = k - (n - 1) * p
        for i in (0, 1):
            wi = (1 - p) if i == 0 else p
            ui = i - p
            si = s_minus + ui
            bi = ui * si - p * (1 - p)
            suffix = str(i)
            count = row["count"]
            assert isinstance(count, F)
            sum_q = allocated(row, f"sum_q{suffix}", n, args.a_raw_denominator)
            sum_e = allocated(row, f"sum_e{suffix}", n, args.a_raw_denominator)
            sum_a = allocated(row, f"sum_a{suffix}", n, args.a_raw_denominator)
            sum_qa = allocated(row, f"sum_q{suffix}a{suffix}", n, args.a_raw_denominator)
            sum_ea = allocated(row, f"sum_e{suffix}a{suffix}", n, args.a_raw_denominator)
            sum_h = 2 * c[geometry] * sum_e - r * sum_q - mu_h[geometry] * count
            sum_ha = (
                2 * c[geometry] * (sum_ea - mu_a * sum_e)
                - r * (sum_qa - mu_a * sum_q)
                - mu_h[geometry] * (sum_a - mu_a * count)
            )
            geometry_cells[(geometry, tau, alpha)] += weight * eq10_state_term(
                wi, ui, bi, sum_h, sum_ha, beta_component, args.z_orbit_multiplicity
            )

    taus = sorted({key[1] for key in geometry_cells})
    alphas = sorted({key[2] for key in geometry_cells})
    expected_cells = {
        (geometry, tau, alpha)
        for geometry in GEOMETRIES
        for tau in taus
        for alpha in alphas
    }
    missing = expected_cells - set(geometry_cells)
    if missing:
        raise ValueError(f"incomplete geometry/tau/alpha rectangle; first missing={sorted(missing)[0]}")

    final_cells = {
        (tau, alpha): (
            geometry_cells[("axis", tau, alpha)]
            + geometry_cells[("tilted", tau, alpha)]
        )
        / 2
        for tau in taus
        for alpha in alphas
    }
    minors = []
    tested_minor_count = 0
    stopped_at_first_nonzero = False
    for tau0, tau1 in combinations(taus, 2):
        for alpha0, alpha1 in combinations(alphas, 2):
            tested_minor_count += 1
            determinant = (
                final_cells[(tau0, alpha0)] * final_cells[(tau1, alpha1)]
                - final_cells[(tau0, alpha1)] * final_cells[(tau1, alpha0)]
            )
            minors.append(
                {
                    "rows": [tau0, tau1],
                    "columns": [alpha0, alpha1],
                    "determinant": interval_record(determinant),
                }
            )
            if args.first_nonzero_only and (determinant.lo > 0 or determinant.hi < 0):
                stopped_at_first_nonzero = True
                break
        if stopped_at_first_nonzero:
            break

    return {
        "schema": "matching-one/p537-siteflip-schur-score/v1",
        "status": "scored" if minors else "matrix_has_fewer_than_two_rows_or_columns",
        "N": n,
        "Delta": str(delta),
        "root_p": interval_record(p),
        "source_allocation": {
            "producer_fields": "raw fixed-source-anchor alpha sums g_0y",
            "scorer_divisor_for_every_a-containing_field": n,
            "raw_a_denominator": args.a_raw_denominator,
            "identity": "a^y=N^-1 g_0y after fixing the source endpoint by translation",
            "no_second_thermal_site_divisor": True,
        },
        "C4_contract": {
            "fixed_NN_z_physical_orbit_multiplicity": args.z_orbit_multiplicity,
            "labels": "tau/alpha are common C4-canonical physical labels; the multiplier restores the fixed-z orbit",
        },
        "global": {
            "M_t": interval_record(mt),
            "Y_t": interval_record(yt),
            "R": interval_record(r),
            "mu_H": {geometry: interval_record(mu_h[geometry]) for geometry in GEOMETRIES},
            "source_components": {
                component: {
                    "jM_component": interval_record(
                        source_packets[("axis", component)]["jM_component"]
                    ),
                    "beta_component": interval_record(beta[component]),
                    "mu_a": {
                        geometry: interval_record(
                            source_packets[(geometry, component)]["mu_a"]
                        )
                        for geometry in GEOMETRIES
                    },
                }
                for component in sorted(beta)
            },
        },
        "matrix": {
            "row_order": taus,
            "column_order": alphas,
            "geometry_cells": {
                geometry: [
                    [interval_record(geometry_cells[(geometry, tau, alpha)]) for alpha in alphas]
                    for tau in taus
                ]
                for geometry in GEOMETRIES
            },
            "P4_Schur": [
                [interval_record(final_cells[(tau, alpha)]) for alpha in alphas] for tau in taus
            ],
        },
        "minors": minors,
        "nonzero_minor_count": sum(item["determinant"]["excludes_zero"] for item in minors),
        "minor_search": {
            "mode": "first_exact_nonzero" if args.first_nonzero_only else "all",
            "tested": tested_minor_count,
            "stopped_at_first_nonzero": stopped_at_first_nonzero,
            "total_possible": math.comb(len(taus), 2) * math.comb(len(alphas), 2),
        },
        "inputs": {
            "aggregates": {"path": str(args.aggregates), "sha256": sha256(args.aggregates)},
            "baseline_axis": {"path": str(args.baseline_axis), "sha256": sha256(args.baseline_axis)},
            "baseline_tilted": {
                "path": str(args.baseline_tilted),
                "sha256": sha256(args.baseline_tilted),
            },
            "baseline_root": {"path": str(args.baseline_root), "sha256": sha256(args.baseline_root)},
        },
    }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--aggregates", type=Path, required=True)
    result.add_argument("--baseline-axis", type=Path, required=True)
    result.add_argument("--baseline-tilted", type=Path, required=True)
    result.add_argument("--baseline-root", type=Path, required=True)
    result.add_argument("--output", type=Path, required=True)
    result.add_argument("--n", type=int, default=25)
    result.add_argument("--delta", default="1152/625")
    result.add_argument(
        "--a-raw-denominator",
        type=int,
        default=1,
        help="raw a-field unit denominator; use 16 when the producer stores g16=16*g",
    )
    result.add_argument(
        "--first-nonzero-only",
        action="store_true",
        help="stop after the first exact nonzero minor; sufficient for the frozen rank-one stop gate",
    )
    result.add_argument(
        "--z-orbit-multiplicity",
        type=int,
        default=4,
        help="physical C4 multiplicity of the fixed thermal-site representative (default: 4)",
    )
    return result


def main() -> None:
    if sys.argv[1:] == ["--self-test"]:
        # p=1/2, source absent (sum_HA=0), only the i=0 Schur cell survives.
        # Its single-z value is -1/4 and the physical C4 orbit is exactly -1.
        term = eq10_state_term(
            Interval.of(F(1, 2)),
            Interval.of(F(-1, 2)),
            Interval.of(F(1, 4)),
            Interval.of(1),
            Interval.of(0),
            Interval.of(2),
            4,
        )
        assert term == Interval.of(-1), term
        single = eq10_state_term(
            Interval.of(F(1, 2)),
            Interval.of(F(-1, 2)),
            Interval.of(F(1, 4)),
            Interval.of(1),
            Interval.of(0),
            Interval.of(2),
            1,
        )
        assert term == 4 * single
        print("self-test passed: absent-source Schur term retained; fixed-z orbit x4")
        return
    args = parser().parse_args()
    payload = score(args)
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
