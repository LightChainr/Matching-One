#!/usr/bin/env python3
"""Local matching-polynomial zeros near the physical critical root.

Issue #113: the global imaginary-RMS cloud of PR #78/#84 is not the finite-size
quantity analogous to a local critical zero. This module reads the committed
exact root CSV and reports only named local diagnostics. It does not fit a
power, does not treat the roots as Fisher/Lee–Yang zeros, and freezes the
metric definitions before any future exact L=6 result is inspected.

Named diagnostics, not fit-and-rescue targets:

    imag_times_L_to_3_over_4      = L^{3/4} |Im z_nearest|
    complex_distance_times_L_to_4 = L^4 |z_nearest - p*|

At the available exact sizes neither diagnostic is stable. The complex-zero
scaling route closes.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CSV = ROOT / "results" / "exact-zero-map-pilot" / "roots.csv"

PHYSICAL_ROOTS = {
    ("axis", 1): 0.5,
    ("axis", 2): 0.54119610014619698439972320536638942006107206337802,
    ("axis", 3): 0.58651145511267563565455897660690173482430062489384,
    ("axis", 4): 0.59067211233102829689590201143951286962111713272216,
    ("diamond", 1): 0.70710678118654752440084436210484903928483593768847,
    ("diamond", 2): 0.60456327785350742069777875145056547767160021550649,
    ("diamond", 3): 0.59425232116856869970538128815606582466846475715026,
}

ALL_REAL = {("axis", 1), ("axis", 2), ("diamond", 1)}


@dataclass(frozen=True)
class ComplexRoot:
    real: float
    imag: float
    is_real: bool
    is_upper_half: bool
    index: int

    @property
    def value(self) -> complex:
        return complex(self.real, self.imag)


@dataclass(frozen=True)
class LocalZeroRow:
    geometry: str
    L: int
    n_roots: int
    n_nonreal: int
    physical_root: float
    nearest_nonreal: Optional[complex]
    imag: Optional[float]
    complex_distance: Optional[float]
    re_in_unit_interval: Optional[bool]
    matching_partner_of_physical: float
    physical_self_matching_gap: float
    local_spacing: float
    nearest_pair_gap: float
    imag_times_L_to_3_over_4: Optional[float]
    complex_distance_times_L_to_4: Optional[float]


def load_roots(path: Path = DEFAULT_CSV) -> dict[tuple[str, int], list[ComplexRoot]]:
    grouped: dict[tuple[str, int], list[ComplexRoot]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "geometry",
            "L",
            "root_index",
            "root_real",
            "root_imaginary",
            "is_real",
            "is_upper_half",
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError("roots CSV missing fields: " + ", ".join(sorted(missing)))
        for raw in reader:
            key = (raw["geometry"], int(raw["L"]))
            grouped.setdefault(key, []).append(
                ComplexRoot(
                    real=float(raw["root_real"]),
                    imag=float(raw["root_imaginary"]),
                    is_real=raw["is_real"] == "True",
                    is_upper_half=raw["is_upper_half"] == "True",
                    index=int(raw["root_index"]),
                )
            )
    if not grouped:
        raise ValueError("roots CSV is empty")
    return grouped


def physical_root(roots: Sequence[ComplexRoot]) -> float:
    physical = [
        root.real
        for root in roots
        if root.is_real and 0.0 < root.real < 1.0 and abs(root.imag) < 1e-18
    ]
    if len(physical) != 1:
        raise ValueError("expected a unique real root in (0,1), got " + str(physical))
    return physical[0]


def nearest_nonreal(roots: Sequence[ComplexRoot], pstar: float) -> Optional[ComplexRoot]:
    nonreal = [root for root in roots if not root.is_real]
    if not nonreal:
        return None
    distance = min(abs(root.value - pstar) for root in nonreal)
    nearest = [root for root in nonreal if abs(abs(root.value - pstar) - distance) <= 1e-15]
    upper = [root for root in nearest if root.is_upper_half]
    chosen = upper[0] if upper else nearest[0]
    return chosen


def nearest_pair_gap(roots: Sequence[ComplexRoot]) -> float:
    values = [root.value for root in roots]
    gap = min(abs(left + right - 1.0) for left in values for right in values)
    return float(gap)


def local_spacing(roots: Sequence[ComplexRoot], pstar: float) -> float:
    others = [abs(root.value - pstar) for root in roots if abs(root.value - pstar) > 1e-12]
    if not others:
        return 0.0
    return float(min(others))


def analyze_polynomial(geometry: str, length: int, roots: Sequence[ComplexRoot]) -> LocalZeroRow:
    pstar = physical_root(roots)
    locked = PHYSICAL_ROOTS[(geometry, length)]
    if abs(pstar - locked) > 1e-12:
        raise AssertionError(
            f"{geometry} L={length}: physical root {pstar} != locked {locked}"
        )
    nearest = nearest_nonreal(roots, pstar)
    imag = None if nearest is None else abs(nearest.imag)
    distance = None if nearest is None else abs(nearest.value - pstar)
    return LocalZeroRow(
        geometry=geometry,
        L=length,
        n_roots=len(roots),
        n_nonreal=sum(0 if root.is_real else 1 for root in roots),
        physical_root=pstar,
        nearest_nonreal=None if nearest is None else nearest.value,
        imag=imag,
        complex_distance=distance,
        re_in_unit_interval=(
            None if nearest is None else 0.0 < nearest.real < 1.0
        ),
        matching_partner_of_physical=1.0 - pstar,
        physical_self_matching_gap=abs(2.0 * pstar - 1.0),
        local_spacing=local_spacing(roots, pstar),
        nearest_pair_gap=nearest_pair_gap(roots),
        imag_times_L_to_3_over_4=(
            None if imag is None else (length ** 0.75) * imag
        ),
        complex_distance_times_L_to_4=(
            None if distance is None else (length ** 4) * distance
        ),
    )


def analyze_catalog(path: Path = DEFAULT_CSV) -> list[LocalZeroRow]:
    grouped = load_roots(path)
    rows = [
        analyze_polynomial(geometry, length, grouped[(geometry, length)])
        for geometry, length in sorted(grouped)
    ]
    return rows


def named_diagnostics(rows: Iterable[LocalZeroRow]) -> dict[str, list[float]]:
    imag = []
    distance = []
    for row in rows:
        if row.imag_times_L_to_3_over_4 is not None:
            imag.append(row.imag_times_L_to_3_over_4)
        if row.complex_distance_times_L_to_4 is not None:
            distance.append(row.complex_distance_times_L_to_4)
    return {
        "imag_times_L_to_3_over_4": imag,
        "complex_distance_times_L_to_4": distance,
    }


def max_min_ratio(values: Sequence[float]) -> Optional[float]:
    if len(values) < 2:
        return None
    lo = min(abs(value) for value in values)
    if lo == 0.0:
        return math.inf
    return max(abs(value) for value in values) / lo


def complex_zero_route_closed(rows: Sequence[LocalZeroRow]) -> bool:
    """Close the route unless a named diagnostic is stable across available sizes.

    Stability here is a descriptive gate, not a fit: every geometry that has
    at least two nonreal sizes must keep both named diagnostics inside a
    factor-of-1.1 envelope. Available exact sizes fail this, and nearest
    nonreal roots are usually outside the physical interval (0,1).
    """

    by_geometry: dict[str, list[LocalZeroRow]] = {}
    for row in rows:
        if row.nearest_nonreal is None:
            continue
        by_geometry.setdefault(row.geometry, []).append(row)
    if not by_geometry:
        return True
    for group in by_geometry.values():
        if len(group) < 2:
            continue
        imag = [row.imag_times_L_to_3_over_4 for row in group]
        dist = [row.complex_distance_times_L_to_4 for row in group]
        imag_ratio = max_min_ratio([value for value in imag if value is not None])
        dist_ratio = max_min_ratio([value for value in dist if value is not None])
        if imag_ratio is not None and imag_ratio <= 1.1 and dist_ratio is not None and dist_ratio <= 1.1:
            return False
    return True


def payload(path: Path = DEFAULT_CSV) -> dict[str, object]:
    rows = analyze_catalog(path)
    diagnostics = named_diagnostics(rows)
    closed = complex_zero_route_closed(rows)
    return {
        "schema": "issue-113 local matching zeros v1",
        "source_csv": str(path.relative_to(ROOT)) if path.is_absolute() else str(path),
        "claim_level": "C1",
        "roots_are": "zeros of the finite matching polynomial",
        "roots_are_not": "Fisher or Lee-Yang zeros",
        "complex_zero_scaling_route": "closed" if closed else "open",
        "named_diagnostics": diagnostics,
        "named_diagnostic_max_min_ratios": {
            name: max_min_ratio(values) for name, values in diagnostics.items()
        },
        "rows": [_row_payload(row) for row in rows],
        "passed": closed and all(
            (row.geometry, row.L) in PHYSICAL_ROOTS for row in rows
        ),
    }


def _row_payload(row: LocalZeroRow) -> dict[str, object]:
    nearest = row.nearest_nonreal
    return {
        "geometry": row.geometry,
        "L": row.L,
        "n_roots": row.n_roots,
        "n_nonreal": row.n_nonreal,
        "physical_root": row.physical_root,
        "nearest_nonreal_real": None if nearest is None else nearest.real,
        "nearest_nonreal_imag": None if nearest is None else nearest.imag,
        "imag": row.imag,
        "complex_distance": row.complex_distance,
        "re_in_unit_interval": row.re_in_unit_interval,
        "matching_partner_of_physical": row.matching_partner_of_physical,
        "physical_self_matching_gap": row.physical_self_matching_gap,
        "local_spacing": row.local_spacing,
        "nearest_pair_gap": row.nearest_pair_gap,
        "imag_times_L_to_3_over_4": row.imag_times_L_to_3_over_4,
        "complex_distance_times_L_to_4": row.complex_distance_times_L_to_4,
        "all_real": nearest is None,
    }


def render_report(data: dict[str, object]) -> str:
    lines = [
        "# Local matching-polynomial zeros near the physical root",
        "",
        "Source: `results/exact-zero-map-pilot/roots.csv` via `scripts/local_matching_zeros.py`.",
        "Claim level: C1 descriptive. These are matching-polynomial zeros, not Fisher zeros.",
        "",
        "Metrics frozen in `predictions/local_matching_zero_metrics_20260829.yaml` before",
        "any future exact L=6 result is inspected. Named diagnostics are not fit targets.",
        "",
        "## Local catalogue",
        "",
        "| geometry | L | physical root | nearest nonreal | `L^{3/4} Im` | `L^4 |z-p*|` | Re in (0,1) |",
        "|---|---:|---:|---|---:|---:|:---:|",
    ]
    for row in data["rows"]:
        if row["all_real"]:
            nearest = "all real"
            imag_diag = "—"
            dist_diag = "—"
            interval = "—"
        else:
            imag_part = row["nearest_nonreal_imag"]
            sign = "+" if imag_part >= 0 else ""
            nearest = "{:.6f}{}{:.6f}i".format(
                row["nearest_nonreal_real"], sign, imag_part
            )
            imag_diag = "{:.3f}".format(row["imag_times_L_to_3_over_4"])
            dist_diag = "{:.1f}".format(row["complex_distance_times_L_to_4"])
            interval = "yes" if row["re_in_unit_interval"] else "no"
        lines.append(
            "| {geometry} | {L} | {physical_root:.12f} | {nearest} | {imag_diag} | {dist_diag} | {interval} |".format(
                geometry=row["geometry"],
                L=row["L"],
                physical_root=row["physical_root"],
                nearest=nearest,
                imag_diag=imag_diag,
                dist_diag=dist_diag,
                interval=interval,
            )
        )
    lines.extend(
        [
            "",
            "## Named-diagnostic stability",
            "",
            "```text",
            "L^{3/4} |Im|  values: "
            + ", ".join("{:.3f}".format(value) for value in data["named_diagnostics"]["imag_times_L_to_3_over_4"]),
            "L^4 |z-p*|   values: "
            + ", ".join("{:.1f}".format(value) for value in data["named_diagnostics"]["complex_distance_times_L_to_4"]),
            "complex-zero scaling route: " + str(data["complex_zero_scaling_route"]),
            "```",
            "",
            "Axis `L^{3/4} Im` jumps from 0.429 (L=3) to 0.730 (L=4). Distance scaled",
            "by `L^4` is 48, 132, 12, 41. Nearest nonreal roots usually have real part",
            "outside `(0,1)`. No stable local power is present at available exact sizes.",
            "",
            "## Boundary",
            "",
            "Do not invent further cloud statistics. Do not score a future L=6 polynomial",
            "against a power fitted here. A later exact size may reopen the route only by",
            "the frozen metrics above, not by a new ad hoc summary.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    data = payload(args.csv)
    status = "CLOSED" if data["complex_zero_scaling_route"] == "closed" else "OPEN"
    print("local matching zeros: complex-zero route " + status)
    for row in data["rows"]:
        extra = "all-real" if row["all_real"] else (
            "L^{{3/4}}Im={imag_times_L_to_3_over_4:.3f} "
            "L^4dist={complex_distance_times_L_to_4:.1f}".format(**row)
        )
        print("{geometry} L={L}: p*={physical_root:.12f} {extra}".format(extra=extra, **row))
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
        print("wrote " + str(args.json))
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(render_report(data), encoding="utf-8")
        print("wrote " + str(args.report))
    return 0 if data["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
