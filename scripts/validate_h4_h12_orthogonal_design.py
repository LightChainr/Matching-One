#!/usr/bin/env python3
"""Validate the exact, data-free Issue 55 orthogonal-design contract."""

from __future__ import annotations

import argparse
from decimal import Decimal, InvalidOperation, localcontext
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence, Tuple

from gaussian_harmonic_arithmetic import harmonic, norm


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "analysis" / "h4_h12_orthogonal_design_manifest.json"
EXPECTED_SCORING_ORDER = [
    "frozen_h4_only_target_with_source_uncertainty_and_future_target_covariance",
    "zero_effect",
    "declared_two_column_h4_h12_model",
    "held_out_third_alias_ratio_only_if_a12_is_resolved",
]
FORBIDDEN_DATA_KEYS = frozenset(
    {
        "observed_mean",
        "observed_se",
        "samples",
        "sample_count",
        "seed",
        "counter_begin",
        "counter_end",
        "covariance",
        "chi2",
        "p_value",
        "h4_score",
        "h12_score",
    }
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _decimal(value: Any, label: str) -> Decimal:
    _require(isinstance(value, str) and value.strip() == value, "%s must be an exact string" % label)
    try:
        parsed = Decimal(value)
    except InvalidOperation as exc:
        raise ValueError("%s is not a decimal" % label) from exc
    _require(parsed.is_finite(), "%s must be finite" % label)
    return parsed


def _fraction(value: Any, label: str) -> Fraction:
    _require(isinstance(value, str) and value.strip() == value, "%s must be an exact string" % label)
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError("%s is not a fraction" % label) from exc


def _walk_forbidden_data(value: Any, path: str = "$") -> None:
    if isinstance(value, Mapping):
        bad = sorted(FORBIDDEN_DATA_KEYS.intersection(value))
        _require(not bad, "%s contains target-data fields: %s" % (path, ",".join(bad)))
        for key, child in value.items():
            _walk_forbidden_data(child, "%s.%s" % (path, key))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _walk_forbidden_data(child, "%s[%d]" % (path, index))


def _pair(value: Any, label: str) -> Tuple[int, int]:
    _require(
        isinstance(value, list)
        and len(value) == 2
        and all(isinstance(component, int) for component in value),
        "%s must be an integer pair" % label,
    )
    return value[0], value[1]


def _rounded_to_declared_digits(actual: Decimal, frozen_text: str, label: str) -> None:
    frozen = _decimal(frozen_text, label)
    quantum = Decimal(1).scaleb(frozen.as_tuple().exponent)
    _require(abs(actual - frozen) <= abs(quantum) / 2, "%s decimal drift" % label)


def _predictions(
    n: int,
    delta_cos4: Fraction,
    amplitude: Decimal,
    amplitude_se: Decimal,
) -> Tuple[Decimal, Decimal]:
    with localcontext() as context:
        context.prec = 60
        exponent = -Decimal(13) / Decimal(8)
        radial = (exponent * Decimal(n).ln()).exp()
        angular = Decimal(delta_cos4.numerator) / Decimal(delta_cos4.denominator)
        return amplitude * angular * radial, amplitude_se * angular * radial


def validate_manifest(manifest: Mapping[str, Any], note_bytes: bytes) -> dict[str, Any]:
    """Return an exact audit summary or raise ValueError on design drift."""

    _walk_forbidden_data(manifest)
    _require(manifest.get("schema") == "matching-one/h4-h12-orthogonal-design/v1", "unknown schema")
    _require(manifest.get("issue") == 55, "wrong issue")
    _require(manifest.get("status") == "design_only_no_target_data", "manifest is not data-free")
    source = manifest.get("source_note", {})
    _require(source.get("path") == "notes/h4-h12-orthogonal-gaussian-design.md", "wrong source note")
    digest = hashlib.sha256(note_bytes).hexdigest()
    _require(source.get("sha256") == digest, "source-note SHA-256 mismatch")

    model = manifest.get("frozen_model", {})
    amplitude = _decimal(model.get("source_amplitude"), "source amplitude")
    amplitude_se = _decimal(model.get("source_amplitude_se"), "source amplitude SE")
    _require(amplitude > 0 and amplitude_se > 0, "source amplitude inputs must be positive")
    _require(model.get("radial_exponent") == "13/8", "radial exponent drift")
    _require(model.get("formula") == "DeltaM=A4*DeltaCos4*N^(-13/8)", "target formula drift")
    _require(model.get("target_mc_uncertainty_included") is False, "target MC uncertainty was fabricated")

    designs = manifest.get("designs")
    _require(isinstance(designs, list) and len(designs) == 2, "exactly two designs required")
    ids = []
    aliases = []
    audited = []
    for index, design in enumerate(designs):
        label = "design[%d]" % index
        design_id = design.get("id")
        _require(isinstance(design_id, str) and design_id, "%s id missing" % label)
        ids.append(design_id)
        n = design.get("N")
        _require(isinstance(n, int) and n > 0, "%s N must be positive" % label)
        first = _pair(design.get("first"), "%s first" % label)
        second = _pair(design.get("second"), "%s second" % label)
        _require(norm(first) == n and norm(second) == n, "%s is not a same-N pair" % label)
        _require(math.gcd(abs(first[0]), abs(first[1])) == 1, "%s first is not primitive" % label)
        _require(math.gcd(abs(second[0]), abs(second[1])) == 1, "%s second is not primitive" % label)
        delta4 = harmonic(first, 1)[0] - harmonic(second, 1)[0]
        delta12 = harmonic(first, 3)[0] - harmonic(second, 3)[0]
        stored4 = _fraction(design.get("delta_cos4"), "%s delta_cos4" % label)
        stored12 = _fraction(design.get("delta_cos12"), "%s delta_cos12" % label)
        stored_alias = _fraction(
            design.get("alias_ratio_cos12_over_cos4"),
            "%s alias ratio" % label,
        )
        _require(delta4 == stored4, "%s signed delta_cos4 drift" % label)
        _require(delta12 == stored12, "%s signed delta_cos12 drift" % label)
        _require(delta4 > 0, "%s signed H4 contrast must be positive" % label)
        _require(delta12 / delta4 == stored_alias, "%s alias ratio drift" % label)
        mean, source_se = _predictions(n, delta4, amplitude, amplitude_se)
        _rounded_to_declared_digits(mean, design.get("h4_only_target_mean"), "%s target mean" % label)
        _rounded_to_declared_digits(
            source_se,
            design.get("source_coefficient_only_se"),
            "%s source-only SE" % label,
        )
        _require(_decimal(design["h4_only_target_mean"], "%s target mean" % label) > 0, "%s target sign drift" % label)
        aliases.append(stored_alias)
        audited.append(
            {
                "id": design_id,
                "N": n,
                "first": list(first),
                "second": list(second),
                "delta_cos4": str(delta4),
                "delta_cos12": str(delta12),
                "alias_ratio": str(stored_alias),
            }
        )
    _require(len(ids) == len(set(ids)), "duplicate design ids")
    _require(aliases[0] < 0 < aliases[1], "alias ratios must have frozen opposite signs")
    _require(all(Fraction(3, 2) < abs(value) < Fraction(2, 1) for value in aliases), "alias magnitudes left design band")

    _require(manifest.get("scoring_order") == EXPECTED_SCORING_ORDER, "scoring order drift")
    protocol = manifest.get("protocol", {})
    for key in (
        "orientation_order_is_signed",
        "pilot_uses_variance_only_to_freeze_sample_count",
        "preserve_batch_covariance",
    ):
        _require(protocol.get(key) is True, "protocol flag must remain true: %s" % key)
    for key in ("add_per_size_harmonic_mixtures", "inspect_target_before_sample_count_freeze"):
        _require(protocol.get(key) is False, "protocol flag must remain false: %s" % key)
    _require(manifest.get("claim_boundary", {}).get("parent_issue") == "remain open", "parent boundary changed")
    return {
        "schema": manifest["schema"],
        "status": "valid_design_only",
        "source_note_sha256": digest,
        "design_count": len(audited),
        "designs": audited,
        "opposite_alias_signs": True,
        "contains_target_data": False,
        "parent_issue": "remain open",
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    args = parser.parse_args(argv)
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    note_path = ROOT / manifest["source_note"]["path"]
    summary = validate_manifest(manifest, note_path.read_bytes())
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
