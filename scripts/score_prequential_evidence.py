#!/usr/bin/env python3
"""Score frozen Gaussian evidence blocks without double counting.

The manifest uses JSON syntax (which is valid YAML) so the scorer remains
stdlib-only.  Scores are accumulated only for one declared primary view per
raw-data group; sensitivities and protocol-history rows remain visible.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


ADDITIVE_ROLE = "primary"


def _cholesky(matrix):
    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix):
        raise ValueError("covariance must be a nonempty square matrix")
    out = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            value = float(matrix[i][j]) - sum(out[i][k] * out[j][k] for k in range(j))
            if i == j:
                if value <= 0.0:
                    raise ValueError("covariance must be positive definite")
                out[i][j] = math.sqrt(value)
            else:
                out[i][j] = value / out[j][j]
    return out


def gaussian_score(residual, covariance):
    """Return chi-square and Gaussian negative log predictive density."""
    if len(residual) != len(covariance):
        raise ValueError("residual and covariance dimensions differ")
    chol = _cholesky(covariance)
    solved = []
    for i, value in enumerate(residual):
        solved.append((float(value) - sum(chol[i][j] * solved[j] for j in range(i))) / chol[i][i])
    chi2 = sum(value * value for value in solved)
    logdet = 2.0 * sum(math.log(chol[i][i]) for i in range(len(chol)))
    nlpd = 0.5 * (chi2 + logdet + len(residual) * math.log(2.0 * math.pi))
    return {"chi_square": chi2, "log_determinant": logdet, "nlpd": nlpd, "dimension": len(residual)}


def validate_manifest(manifest):
    seen_ids = set()
    primary_groups = {}
    for block in manifest["blocks"]:
        block_id = block["id"]
        if block_id in seen_ids:
            raise ValueError(f"duplicate block id: {block_id}")
        seen_ids.add(block_id)
        if block.get("role") == ADDITIVE_ROLE and block.get("status") == "SCORED":
            group = block["raw_data_group"]
            if group in primary_groups:
                raise ValueError(
                    f"raw-data group {group} has multiple additive primary views: "
                    f"{primary_groups[group]} and {block_id}"
                )
            primary_groups[group] = block_id
        channel = block.get("channel", {})
        if channel.get("source") != channel.get("target") and not channel.get("exact_map"):
            if block.get("status") == "SCORED":
                raise ValueError(f"scored block {block_id} has an unregistered channel mismatch")


def score_manifest(manifest):
    validate_manifest(manifest)
    rows = []
    cumulative = {}
    coverage = {}
    for block in manifest["blocks"]:
        row = {
            "id": block["id"], "raw_data_group": block["raw_data_group"],
            "role": block["role"], "status": block["status"],
            "channel": block.get("channel"), "scores": {},
        }
        if block["status"] == "SCORED":
            observation = block["observation"]
            for model, prediction in block["models"].items():
                coverage.setdefault(model, {"predicted": [], "additive": []})["predicted"].append(block["id"])
                if prediction.get("status") == "NO_PREDICTION":
                    continue
                residual = [x - m for x, m in zip(observation, prediction["mean"])]
                score = gaussian_score(residual, prediction["covariance"])
                row["scores"][model] = score
                if block["role"] == ADDITIVE_ROLE:
                    coverage[model]["additive"].append(block["id"])
                    total = cumulative.setdefault(model, {"chi_square": 0.0, "nlpd": 0.0, "dimensions": 0, "blocks": []})
                    total["chi_square"] += score["chi_square"]
                    total["nlpd"] += score["nlpd"]
                    total["dimensions"] += score["dimension"]
                    total["blocks"].append(block["id"])
        rows.append(row)

    pairwise = []
    models = sorted(cumulative)
    score_by_id = {row["id"]: row["scores"] for row in rows}
    for i, left in enumerate(models):
        for right in models[i + 1:]:
            intersection = sorted(set(coverage[left]["additive"]) & set(coverage[right]["additive"]))
            if not intersection:
                continue
            left_score = sum(score_by_id[b][left]["nlpd"] for b in intersection)
            right_score = sum(score_by_id[b][right]["nlpd"] for b in intersection)
            pairwise.append({
                "left": left, "right": right, "intersection": intersection,
                "delta_nlpd_left_minus_right": left_score - right_score,
                "preferred": left if left_score < right_score else right,
            })
    return {
        "schema_version": 1,
        "manifest_version": manifest["schema_version"],
        "blocks": rows,
        "cumulative_primary_only": cumulative,
        "coverage": coverage,
        "pairwise_primary_intersections": pairwise,
        "governance": manifest["governance"],
    }


def render_markdown(result):
    lines = ["# Prequential evidence ledger", "", "Only `primary` rows enter cumulative evidence.", "",
             "| Block | Role | Status | Model | chi2 | NLPD |", "|---|---|---|---|---:|---:|"]
    for block in result["blocks"]:
        if not block["scores"]:
            lines.append(f"| {block['id']} | {block['role']} | {block['status']} | — | — | — |")
        for model, score in sorted(block["scores"].items()):
            lines.append(f"| {block['id']} | {block['role']} | {block['status']} | {model} | {score['chi_square']:.6g} | {score['nlpd']:.6f} |")
    lines += ["", "## Pairwise comparisons on matched primary endpoints", "",
              "Negative delta favors the left model.", "", "| Left | Right | Blocks | Delta NLPD | Preferred |",
              "|---|---|---|---:|---|"]
    for row in result["pairwise_primary_intersections"]:
        lines.append(f"| {row['left']} | {row['right']} | {', '.join(row['intersection'])} | {row['delta_nlpd_left_minus_right']:.6f} | {row['preferred']} |")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--markdown-output", type=Path)
    args = parser.parse_args()
    result = score_manifest(json.loads(args.manifest.read_text()))
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.json_output:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(rendered)
    else:
        print(rendered, end="")
    if args.markdown_output:
        args.markdown_output.parent.mkdir(parents=True, exist_ok=True)
        args.markdown_output.write_text(render_markdown(result))


if __name__ == "__main__":
    main()
