#!/usr/bin/env python3
"""Typed generic-Q lift descriptors and exact first-tangent transport."""

from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Iterable, Mapping


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PHASE_A = ROOT / "results/q-lift-covariance/latest.json"
DEFAULT_OUTPUT = ROOT / "analysis/generic_q_lift_transport_manifest.json"
ENDPOINT = "A_top_Q1=P_2D-P_0D"
SECTOR_ORDER = ("2D", "1D", "0D")


PATHS: dict[str, dict[str, Any]] = {
    "fixed_v_1": {
        "parameter": "Q",
        "Q": "Q",
        "v": "1",
        "endpoint": {"Q": "1", "v": "1"},
        "dv_dQ_at_1": "0",
    },
    "critical_square_bond_v_sqrt_Q": {
        "parameter": "Q",
        "Q": "Q",
        "v": "sqrt(Q)",
        "endpoint": {"Q": "1", "v": "1"},
        "dv_dQ_at_1": "1/2",
        "critical_constraint": "W_2D(Q,v(Q))=Q W_0D(Q,v(Q))",
    },
}


LIFTS: dict[str, dict[str, Any]] = {
    "L_hom": {
        "sector_weights_in_Q": {"2D": "1", "1D": "0", "0D": "-1"},
        "sector_weight_derivative_at_Q1": {"2D": "0", "1D": "0", "0D": "0"},
        "explicit_insertion_Q_dependence": "none; the ambient homology-rank weights are Q-independent",
        "projector_convention": "unweighted ambient-homology projector P_2D-P_0D",
        "counterterm_convention": "raw homology section; CP-horizontal transport subtracts the endpoint 0D sector",
    },
    "L_CP": {
        "sector_weights_in_Q": {"2D": "1", "1D": "0", "0D": "-Q"},
        "sector_weight_derivative_at_Q1": {"2D": "0", "1D": "0", "0D": "-1"},
        "explicit_insertion_Q_dependence": "the 0D projector carries the explicit coefficient -Q",
        "projector_convention": "critical-polynomial / periodic-TL projector P_2D-Q P_0D",
        "counterterm_convention": "CP-horizontal reference section on the exact critical manifold",
    },
}

NORMALIZATIONS: dict[str, dict[str, str]] = {
    "normalized_probability": {
        "numerator": "sum_r w_r(Q) W_r(Q,v)",
        "denominator": "Z(Q,v)=W_0D+W_1D+W_2D",
        "endpoint_sector_basis": "pi_rD=W_rD/Z",
    },
    "restricted_state_sum": {
        "numerator": "sum_r w_r(Q) W_r(Q,v)",
        "denominator": "1",
        "endpoint_sector_basis": "W_rD",
    },
}


@dataclass(frozen=True)
class GenericQLiftDescriptor:
    endpoint_observable_id: str
    lift_id: str
    normalization: str
    path_id: str
    field_normalization_convention: str = "no extra field rescaling; sector normalization only"

    def __post_init__(self) -> None:
        if self.endpoint_observable_id != ENDPOINT:
            raise ValueError("unknown Q=1 endpoint observable")
        if self.lift_id not in LIFTS:
            raise ValueError("unknown generic-Q lift")
        if self.normalization not in NORMALIZATIONS:
            raise ValueError("unknown normalization")
        if self.path_id not in PATHS:
            raise ValueError("unknown (Q,v) path")

    @property
    def descriptor_id(self) -> str:
        return "/".join((self.endpoint_observable_id, self.lift_id, self.normalization, self.path_id))

    def to_dict(self) -> dict[str, Any]:
        lift = LIFTS[self.lift_id]
        return {
            "schema": "matching-one/generic-q-lift-descriptor/v1",
            "descriptor_id": self.descriptor_id,
            "endpoint_observable_id": self.endpoint_observable_id,
            "lift_id": self.lift_id,
            "sector_order": list(SECTOR_ORDER),
            "sector_weights_in_Q": copy.deepcopy(lift["sector_weights_in_Q"]),
            "normalization": {
                "id": self.normalization,
                **copy.deepcopy(NORMALIZATIONS[self.normalization]),
            },
            "Q_v_path": {"id": self.path_id, **copy.deepcopy(PATHS[self.path_id])},
            "explicit_insertion_Q_dependence": lift["explicit_insertion_Q_dependence"],
            "insertion_weight_derivative_at_Q1": copy.deepcopy(lift["sector_weight_derivative_at_Q1"]),
            "projector_convention": lift["projector_convention"],
            "counterterm_convention": lift["counterterm_convention"],
            "field_normalization_convention": self.field_normalization_convention,
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "GenericQLiftDescriptor":
        if payload.get("schema") != "matching-one/generic-q-lift-descriptor/v1":
            raise ValueError("unknown generic-Q lift descriptor schema")
        normalization = payload.get("normalization")
        path = payload.get("Q_v_path")
        if not isinstance(normalization, Mapping) or not isinstance(path, Mapping):
            raise ValueError("descriptor lacks normalization or (Q,v) path")
        field_convention = payload.get("field_normalization_convention")
        if not isinstance(field_convention, str) or not field_convention:
            raise ValueError("descriptor lacks field normalization/counterterm convention")
        value = cls(
            str(payload.get("endpoint_observable_id")),
            str(payload.get("lift_id")),
            str(normalization.get("id")),
            str(path.get("id")),
            field_convention,
        )
        canonical = value.to_dict()
        for key in (
            "descriptor_id", "sector_order", "sector_weights_in_Q",
            "normalization", "Q_v_path", "explicit_insertion_Q_dependence",
            "insertion_weight_derivative_at_Q1", "projector_convention",
            "counterterm_convention",
        ):
            if payload.get(key) != canonical[key]:
                raise ValueError(f"generic-Q descriptor drift: {key}")
        return value


def descriptor(
    lift_id: str,
    path_id: str = "critical_square_bond_v_sqrt_Q",
    normalization: str = "normalized_probability",
) -> GenericQLiftDescriptor:
    return GenericQLiftDescriptor(ENDPOINT, lift_id, normalization, path_id)


def raw_tangents_directly_comparable(
    left: GenericQLiftDescriptor, right: GenericQLiftDescriptor
) -> bool:
    """Raw tangents are directly comparable only for identical descriptors."""

    return left == right


def require_raw_tangent_identity(
    left: GenericQLiftDescriptor, right: GenericQLiftDescriptor
) -> None:
    if not raw_tangents_directly_comparable(left, right):
        raise ValueError(
            "raw Q tangents are not directly comparable: endpoint equality does not "
            "erase lift/path/normalization semantics; apply a registered transport"
        )


def _weight_derivatives(lift_id: str) -> tuple[Fraction, Fraction, Fraction]:
    values = LIFTS[lift_id]["sector_weight_derivative_at_Q1"]
    return tuple(Fraction(values[sector]) for sector in SECTOR_ORDER)  # type: ignore[return-value]


def transport_tangent(
    source: GenericQLiftDescriptor,
    target: GenericQLiftDescriptor,
    *,
    endpoint_sectors: Mapping[str, Fraction | int | str] | None = None,
    tangent: Fraction | int | str | None = None,
) -> dict[str, Any]:
    """Return the exact lift transition at fixed path and normalization.

    The registered map is induced solely by the explicit derivative of the
    sector weights.  At Q=1 both lifts have identical endpoint weights, so all
    measure, path and denominator derivatives cancel from their difference.
    """

    if source.endpoint_observable_id != target.endpoint_observable_id:
        raise ValueError("endpoint observables differ")
    if source.path_id != target.path_id:
        raise ValueError("no transport across different (Q,v) paths")
    if source.normalization != target.normalization:
        raise ValueError("no transport across different normalizations")
    if source.field_normalization_convention != target.field_normalization_convention:
        raise ValueError("field normalization/counterterm convention differs")
    source_derivative = _weight_derivatives(source.lift_id)
    target_derivative = _weight_derivatives(target.lift_id)
    shift_coefficients = tuple(
        target_derivative[i] - source_derivative[i] for i in range(3)
    )
    basis_prefix = "pi_" if source.normalization == "normalized_probability" else "W_"
    basis = tuple(f"{basis_prefix}{sector}" for sector in SECTOR_ORDER)
    expression_terms = []
    for coefficient, name in zip(shift_coefficients, basis):
        if coefficient:
            expression_terms.append(f"{coefficient}*{name}")
    shift_expression = "0" if not expression_terms else "+".join(expression_terms).replace("-1*", "-").replace("1*", "")
    result: dict[str, Any] = {
        "schema": "matching-one/generic-q-tangent-transport/v1",
        "source_descriptor_id": source.descriptor_id,
        "target_descriptor_id": target.descriptor_id,
        "endpoint_weights_equal": True,
        "path_held_fixed": source.path_id,
        "normalization_held_fixed": source.normalization,
        "shift_basis_order": list(basis),
        "shift_coefficients": [str(value) for value in shift_coefficients],
        "tangent_shift_target_minus_source": shift_expression,
        "derivation": "difference of explicit insertion-weight derivatives at Q=1; measure/path/normalization derivatives cancel",
        "exact": True,
    }
    if endpoint_sectors is not None:
        values = tuple(Fraction(endpoint_sectors[name]) for name in basis)
        numeric_shift = sum(
            (coefficient * value for coefficient, value in zip(shift_coefficients, values)),
            Fraction(0),
        )
        result["numeric_shift"] = str(numeric_shift)
        if tangent is not None:
            result["source_tangent"] = str(Fraction(tangent))
            result["transported_tangent"] = str(Fraction(tangent) + numeric_shift)
    elif tangent is not None:
        raise ValueError("a tangent value requires endpoint sector coordinates")
    return result


def to_cp_horizontal(
    source: GenericQLiftDescriptor,
    *,
    endpoint_sectors: Mapping[str, Fraction | int | str],
    tangent: Fraction | int | str,
) -> dict[str, Any]:
    target = descriptor("L_CP", source.path_id, source.normalization)
    result = transport_tangent(
        source, target, endpoint_sectors=endpoint_sectors, tangent=tangent
    )
    result["connection"] = "CP-horizontal"
    result["comparison_scope"] = "transported tangents from #258/#262/#263/#275 may be compared only after their complete descriptors pass"
    return result


def _fraction(value: Mapping[str, Any]) -> Fraction:
    return Fraction(int(value["numerator"]), int(value["denominator"]))


def build_manifest(phase_a: Mapping[str, Any]) -> dict[str, Any]:
    if phase_a.get("schema") != "matching-one/q-lift-covariance-oracle/v1":
        raise ValueError("wrong Phase-A oracle")
    descriptors = [
        descriptor(lift, path, normalization).to_dict()
        for normalization in NORMALIZATIONS
        for path in PATHS
        for lift in LIFTS
    ]
    checks = []
    for row in phase_a["finite_tori"]:
        pi0 = _fraction(row["Q1_v1_sector_counts"]["pi_0D"])
        endpoint = {"pi_2D": pi0, "pi_1D": 1 - 2 * pi0, "pi_0D": pi0}
        for path_id, path in row["paths"].items():
            source = descriptor("L_hom", path_id)
            target = descriptor("L_CP", path_id)
            d_h = _fraction(path["normalized_dh"])
            d_c = _fraction(path["normalized_dc"])
            transport = transport_tangent(
                source, target, endpoint_sectors=endpoint, tangent=d_h
            )
            passed = (
                transport["numeric_shift"] == str(-pi0)
                and transport["transported_tangent"] == str(d_c)
            )
            checks.append({
                "L": row["L"],
                "path_id": path_id,
                "pi_0D": str(pi0),
                "raw_d_L_hom": str(d_h),
                "raw_d_L_CP": str(d_c),
                "transport": transport,
                "passed": passed,
            })
    return {
        "schema": "matching-one/generic-q-lift-transport-manifest/v1",
        "issue": 333,
        "phase": "B",
        "endpoint_observable_id": ENDPOINT,
        "registered_descriptors": descriptors,
        "registered_connection": {
            "id": "CP-horizontal",
            "horizontal_lift": "L_CP",
            "normalized_L_hom_to_L_CP_tangent_shift": "-pi_0D",
            "unnormalized_L_hom_to_L_CP_tangent_shift": "-W_0D",
            "claim": "same endpoint plus different lift forbids direct raw-tangent pooling; compare only after exact CP-horizontal transport",
        },
        "phase_A_transport_checks": checks,
        "all_checks_passed": all(row["passed"] for row in checks),
        "affected_programs": [258, 262, 263, 275],
        "scientific_boundary": (
            "This registers one exact comparison connection; it neither declares CP-horizontal "
            "uniquely physical nor supplies missing representation-projector or field counterterms."
        ),
    }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase-a", type=Path, default=DEFAULT_PHASE_A)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(list(argv) if argv is not None else None)
    phase_a = json.loads(args.phase_a.read_text(encoding="utf-8"))
    result = build_manifest(phase_a)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, allow_nan=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(args.output)
    return 0 if result["all_checks_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
