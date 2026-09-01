#!/usr/bin/env python3
"""Mixed-precision source response on nested old norm-4 production subsets."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import subprocess
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import scipy

import analyze_norm4_source_thermal as old
from analyze_norm4_source_two_phase import covariance_component
from norm4_source_two_phase_core import baseline, _response, _topology_response
from norm4_source_two_phase_inputs import load_complement

ROOT, NS = old.ROOT, old.NS
OUTPUT = ROOT / "results/norm4-source-endpoint-1m"
ENDPOINTS = (260, 340)
EXTRA_FIELDS = ("root_comoving_rank1_fugacity", "rootdot_first_activation_fugacity",
                "rootdot_second_completion_fugacity", "root_comoving_E_H4_fugacity",
                "rank1_common_q", "E_H4_common_q", "Udot_common_q",
                "rank1_common_E", "E_H4_common_E", "Udot_common_E",
                "det_common_q", "det_common_E", "det_common_q_E")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_profile(path, n, expected_samples, run):
    """Read equal-size paired batches; source coordinates return in s/N units."""
    values = np.zeros((100, 2, n + 1, 5), dtype=np.float64)
    seen = set()
    with path.open() as handle:
        for row in csv.DictReader(handle):
            b, k = int(row["batch"]), int(row["k"])
            g = ("first", "second").index(row["orientation"])
            key = (b, g, k)
            if (key in seen or int(row["n"]) != n or not 0 <= b < 100
                    or not 0 <= k <= n or int(row["samples"]) != expected_samples
                    or [int(row["a"]), int(row["b"])] != run[("first", "second")[g]]):
                raise ValueError(f"{path}: inconsistent or duplicate profile row {key}")
            seen.add(key)
            values[b, g, k] = [int(row[f]) for f in
                               ("sum_q", "sum_e", "sum_s", "sum_qs", "sum_es")]
    if len(seen) != 100 * 2 * (n + 1):
        raise ValueError(f"{path}: incomplete all-K paired profile")
    values[..., 2:] /= n
    return values


def extra_readouts(direction, n):
    rows = [direction[name] for name in ("first", "second")]
    d = sum(row["q_p"] for row in rows) / 2
    jq = sum(row["Jq"] for row in rows) / 2
    je = sum(row["JE"] for row in rows) / 2
    ep = sum(row["E_p"] for row in rows) / 2
    result = dict(zip(EXTRA_FIELDS[:3], map(float, (
        -n * (je - ep * jq / d), -n * (jq - je) / (2 * d),
        -n * (jq + je) / (2 * d)))))
    delta = float(Fraction("1152/845" if n in (65, 130, 260) else "2304/1445"))
    anchor = {key: [row[field] for row in rows] for key, field in
              (("m", "q"), ("e", "E"), ("m_p", "q_p"), ("e_p", "E_p"))}
    anchor.update(D=d, B=(rows[0]["E_p"] - rows[1]["E_p"]) / delta,
                  H=(rows[0]["E_pp"] - rows[1]["E_pp"]) / delta,
                  T=(rows[0]["q_pp"] + rows[1]["q_pp"]) / 2)
    source = _response(anchor, n, delta, *[[row[field] for row in rows]
                       for field in ("Jq", "JE", "Jq_p", "JE_p")])

    def column(response, scale):
        qj = np.asarray(response["Jq_density"])
        ej = np.asarray(response["JE_density"])
        moving_e = ej - np.asarray(anchor["e_p"]) * qj.mean() / d
        return scale * np.array([-moving_e.mean(), (moving_e[0] - moving_e[1]) / delta,
                                 response["density"]])

    s_column = column(source, n)
    result["root_comoving_E_H4_fugacity"] = float(s_column[1])
    columns = [s_column]
    for name, indices in (("q", (3, 4)), ("E", (5, 6))):
        coefficients = np.zeros(7)
        coefficients[list(indices)] = 1
        # These columns couple to exp(lambda*q) or exp(eta*E), WITHOUT factor N.
        control = column(_topology_response(coefficients, anchor, n, delta), 1)
        columns.append(control)
        for key, value in zip((f"rank1_common_{name}", f"E_H4_common_{name}",
                               f"Udot_common_{name}"), control):
            result[key] = float(value)
        result[f"det_common_{name}"] = float(s_column[0] * control[2] - control[0] * s_column[2])
    result["det_common_q_E"] = float(np.linalg.det(np.column_stack(columns)))
    return result


def extras_at_saved_root(sums, samples, n, p0):
    direction = {name: old.direction_values(old.binomial_moments(sums[g], samples, p0, n))
                 for g, name in enumerate(("first", "second"))}
    return extra_readouts(direction, n)


def vectorize(points, anchors):
    """Raw U/source remains primary; generator drift uses independent U anchors."""
    values = old.vectorize(points)
    for n in NS:
        values.update({f"N{n}.{field}": float(points[n][field]) for field in EXTRA_FIELDS})
        values[f"N{n}.U_anchor"] = float(anchors[n]["U"])
        values[f"N{n}.p0_anchor"] = float(anchors[n]["p0"])
    for start, middle, end in old.LINEAGES:
        values[f"drift_shape_within_source_supplement.{start}"] = values[f"drift_shape.{start}"]
        values[f"drift_shape.{start}"] = anchors[middle]["U"] - anchors[start]["U"]
        for model, coefficients in (("q2", (1, -3, 2)), ("Jordan", (1, -2, 1))):
            values[f"{model}.{start}.U_anchor"] = sum(
                c * anchors[n]["U"] for c, n in zip(coefficients, (start, middle, end)))
    for model, factor in (("q2", .5), ("Jordan", 1.0)):
        key = f"{model}.generator_drift_determinant"
        values[f"{model}.generator_drift_within_source_supplement"] = values[key]
        values[key] = factor * (
            values[f"{model}.65.Udot_fugacity"] * values["drift_shape.85"]
            - values[f"{model}.85.Udot_fugacity"] * values["drift_shape.65"])
    return {key: float(value) for key, value in values.items()}


def estimate_rows(labels, values, covariance):
    errors = np.sqrt(np.maximum(0.0, np.diag(covariance)))
    return {label: {"value": float(value), "se": float(error),
                    "z": float(value / error) if error else None}
            for label, value, error in zip(labels, values, errors)}


def main():
    started = time.perf_counter()
    destination = OUTPUT / "latest.json"
    if destination.exists():
        raise ValueError("saved endpoint-1M analysis exists; do not overwrite or repeat")
    previous_path = old.OUTPUT / "latest.json"
    previous = json.loads(previous_path.read_text())
    previous_index = {label: i for i, label in enumerate(previous["labels"])}
    manifest = json.loads(old.MANIFEST.read_text())
    hypotheses = json.loads(old.HYPOTHESES.read_text())
    runs = {run["N"]: run for run in manifest["runs"]}
    bracket = hypotheses["root_bracket"]
    deltas = {n: float(Fraction(runs[n]["delta_cos4"])) for n in NS}
    profiles, original, increments, complements, anchors = {}, {}, {}, {}, {}
    old_points, points, diagnostics, groups, inputs = {}, {}, {}, {}, []
    per_batch = {n: 10000 if n in ENDPOINTS else 1000 for n in NS}

    for n in NS:
        old_path = old.OUTPUT / "raw" / f"n{n}.csv"
        original[n] = load_profile(old_path, n, 1000, runs[n])
        inputs.append({"N": n, "role": "original_100k_marks", "path": str(old_path.relative_to(ROOT)),
                       "sha256": digest(old_path)})
        profiles[n] = original[n]
        complements[n] = load_complement(n, original[n])
        if n in ENDPOINTS:
            path = OUTPUT / "increment" / "raw" / f"n{n}.csv"
            increments[n] = load_profile(path, n, 9000, runs[n])
            profiles[n] = original[n] + increments[n]
            inputs.append({"N": n, "role": "disjoint_900k_increment_marks", "path": str(path.relative_to(ROOT)),
                           "sha256": digest(path)})
            complement = complements[n]
            complement["sums"][0] -= increments[n][..., :2].sum(axis=0).astype(np.int64)
            complement["counts"][0] -= 900000
            complement["provenance"].update({
                "removedcount": 1000000, "complementcount": int(complement["counts"].sum()),
                "batchcounts": complement["counts"].tolist(),
                "removed_counter_interval": [8200000000, 8201000000],
                "additional_removed_counter_interval": [8200100000, 8201000000],
                "additional_marked_path": str(path.relative_to(ROOT)), "additional_marked_sha256": digest(path),
                "operation": "Original full1B minus old100k minus disjoint900k; never reuse the overlapping old two-phase complement"})
        complement = complements[n]
        anchors[n] = baseline(complement["sums"].sum(axis=0), int(complement["counts"].sum()),
                              n, deltas[n], bracket)
        old_points[n] = {field: previous["estimates"][f"N{n}.{field}"]["value"] for field in old.FIELDS}
        old_points[n].update(extra_readouts(previous["by_N"][str(n)]["direction"], n))
        if n in ENDPOINTS:
            point, diagnostics[n] = old.at_root(profiles[n].sum(axis=0), 1000000, n, deltas[n], bracket)
            point.update(extra_readouts(diagnostics[n]["direction"], n))
            points[n] = point
        else:
            points[n], diagnostics[n] = dict(old_points[n]), previous["by_N"][str(n)]
        groups.setdefault(runs[n]["dependency_group"], []).append(n)

    totals = {n: profiles[n].sum(axis=0) for n in NS}
    old_totals = {n: original[n].sum(axis=0) for n in NS}
    complement_totals = {n: complements[n]["sums"].sum(axis=0) for n in NS}
    current_map = vectorize(points, anchors)
    current_labels = list(current_map)
    width = len(current_labels)
    labels = current_labels + ["old100k." + label for label in current_labels]

    def packed(current, earlier, current_anchors):
        return np.array(list(vectorize(current, current_anchors).values())
                        + list(vectorize(earlier, current_anchors).values()), dtype=float)

    central = packed(points, old_points, anchors)
    covariance = np.zeros((len(labels), len(labels)))
    contributions = {}
    for group, sizes in groups.items():
        vectors = []
        saved = previous["dependency_groups"][group]["delete_one_vectors"]
        for batch in range(100):
            changed, earlier = dict(points), dict(old_points)
            for n in sizes:
                prior = {field: saved[batch][previous_index[f"N{n}.{field}"]] for field in old.FIELDS}
                prior.update(extras_at_saved_root(old_totals[n] - original[n][batch], 99000, n, prior["p0"]))
                earlier[n] = prior
                if n in ENDPOINTS:
                    state, diagnostic = old.at_root(totals[n] - profiles[n][batch], 990000,
                                                    n, deltas[n], bracket)
                    state.update(extra_readouts(diagnostic["direction"], n))
                    changed[n] = state
                else:
                    changed[n] = dict(prior)
            vectors.append(packed(changed, earlier, anchors))
        vectors = np.asarray(vectors)
        counts = np.repeat(per_batch[sizes[0]], 100)
        component = covariance_component(central, vectors, counts)
        covariance += component
        contributions["source:" + group] = {
            "stage": "marked_source", "Ns": sizes, "batch_counts": counts.tolist(),
            "old100k_batch_counts": [1000] * 100, "delete_one_batch_ids": list(range(100)),
            "delete_one_vectors": vectors.tolist(), "covariance": component.tolist(),
            "operation": "Omit one equal-weight union batch from current marks and its same old1000 subset from the saved old100k estimator; cyclic old U/source vectors reused; all added observables use paired leave-one-out roots"}
        print(f"Completed paired marked-source covariance: {sizes}", flush=True)

    for group, sizes in groups.items():
        counts = complements[sizes[0]]["counts"]
        if not all(np.array_equal(counts, complements[n]["counts"]) for n in sizes):
            raise ValueError("complement batch sizes differ within a shared-counter group")
        vectors = []
        for batch in range(100):
            changed_anchors = dict(anchors)
            for n in sizes:
                complement = complements[n]
                changed_anchors[n] = baseline(complement_totals[n] - complement["sums"][batch],
                    int(complement["counts"].sum() - complement["counts"][batch]), n, deltas[n], bracket)
            vectors.append(packed(points, old_points, changed_anchors))
        vectors = np.asarray(vectors)
        component = covariance_component(central, vectors, counts)
        covariance += component
        contributions["complement:" + group] = {
            "stage": "unmarked_disjoint_complement", "Ns": sizes, "batch_counts": counts.tolist(),
            "delete_one_batch_ids": list(range(100)), "delete_one_vectors": vectors.tolist(),
            "covariance": component.tolist(),
            "operation": "Actual remaining original production batch omitted; anchor U/root and primary generator drift updated, all raw source estimates unchanged; endpoints exclude the entire1M marked union"}
        print(f"Completed independent anchor covariance: {sizes}", flush=True)

    estimates = estimate_rows(labels, central, covariance)
    difference = central[:width] - central[width:]
    difference_covariance = (covariance[:width, :width] + covariance[width:, width:]
                             - covariance[:width, width:] - covariance[width:, :width])
    difference_estimates = estimate_rows(current_labels, difference, difference_covariance)
    models = {}
    for model in ("q2", "Jordan"):
        models[model] = {
            "source_rigidity": old.joint_zero(labels, central, covariance,
                [f"{model}.{n}.Udot_fugacity" for n in (65, 85)]),
            "unperturbed_complement": old.joint_zero(labels, central, covariance,
                [f"{model}.{n}.U_anchor" for n in (65, 85)]),
            "generator_drift_independent_anchor_primary": estimates[f"{model}.generator_drift_determinant"],
            "generator_drift_within_source_supplement": estimates[f"{model}.generator_drift_within_source_supplement"]}

    source_paths = [Path(__file__), Path(old.__file__), ROOT / "scripts/analyze_p40_source_thermal.py",
        ROOT / "scripts/analyze_norm4_source_two_phase.py", ROOT / "scripts/norm4_source_two_phase_core.py",
        ROOT / "scripts/norm4_source_two_phase_inputs.py", old.MANIFEST, old.HYPOTHESES,
        ROOT / "analysis/norm4_source_endpoint_1m_contract.json",
        ROOT / "notes/norm4-root-comoving-rank1-source.md",
        ROOT / "notes/norm4-common-source-response-determinants.md"]
    receipts = [path for path in (OUTPUT / "increment/run.json", OUTPUT / "run.json") if path.exists()]
    result = {
        "schema": "matching-one.norm4-source-endpoint-1m.v1",
        "status": "computed_nested_old_counter_source_extension_mixed_precision",
        "execution_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "labels": labels, "estimates": estimates, "covariance": covariance.tolist(),
        "covariance_contributions": contributions,
        "model_extension_diagnostics": models,
        "by_N": {n: {"source": diagnostics[n], "points": points[n], "anchor": anchors[n],
                     "marked_samples": 100 * per_batch[n]} for n in NS},
        "source_hypotheses": hypotheses,
        "source_coordinate": "s=CB+CW coupled by common microscopic t; v=N*Udot_density; raw original-U source estimator is primary",
        "primary_generator_drift": "R_kappa(v)_65*(U170_anchor-U85_anchor)-R_kappa(v)_85*(U130_anchor-U65_anchor), kappa=.5 or1; every anchor uses full production minus all current marked counters and propagates its own remaining-block covariance",
        "exploratory_clock_null": {
            "rank1": "-N*(barJE-barEp*barJq/D), evaluated at the source estimator's pooled root; not the earlier h*q-compensated C",
            "first_activation_rootdot": "-N*(barJq-barJE)/(2D)",
            "second_completion_rootdot": "-N*(barJq+barJE)/(2D)",
            "boundary": "A nonzero comoving rank1 response rejects pure common thermal reparameterization at that finite root, not a claim of energy-operator identity; the two rootdot allocations are static cumulative-event responses, not marked K1/K2 path attribution"},
        "exploratory_common_source_determinants": {
            "rows": ["pooled root-comoving rank1", "root-comoving P4[E]", "original Udot"],
            "columns": ["bulk s=CB+CW", "q coupled by exp(lambda*q)", "E coupled by exp(eta*E)"],
            "two_by_two": "det_common_q/E uses rows rank1,U; each zero is necessary for its named common single-source plus thermal-clock model",
            "three_by_three": "det_common_q_E uses all three rows; zero is necessary for s mapping to a_g+b*K+c*q+d*E with b,c,d common across the two directions and fixed under thermal differentiation, but allowed to depend on N",
            "boundary": "Exploratory comparisons fixed before aggregation, not selected by p. Nonzero rejects the declared finite first-order common-source map in these readouts, not geometry-dependent source coefficients, a universal three-state theorem or a continuum energy-field identity"},
        "paired_old_vs_new": {
            "labels": current_labels, "difference_direction": "current minus old100k",
            "estimates": difference_estimates, "covariance": difference_covariance.tolist(),
            "joint_vector_order": "Current labels followed by old100k-prefixed copies; every source and anchor leave-one-out vector is saved in that same order",
            "boundary": "Old100k is nested inside the endpoint1M union, not an independent replicate. Old and current generator drift comparisons use the same newly disjoint anchor, not the overlapping previously published two-phase anchor"},
        "sampling": {
            "marked_samples_by_N": {n: 100 * per_batch[n] for n in NS},
            "analysis_batch_counts_by_N": {n: [per_batch[n]] * 100 for n in NS},
            "endpoint_old_counter_interval": [8200000000, 8200100000],
            "endpoint_added_counter_interval": [8200100000, 8201000000],
            "endpoint_batch_union": "old1000 counters in batch b plus disjoint increment9000 counters in batch b, yielding100 equal10000-permutation paired-direction batches",
            "new_random_samples": 0, "newly_marked_old_permutations_this_increment": 1800000,
            "precision_boundary": "Cyclic sizes retain100k source marks; endpoints have1M source marks, not1.9B/1B source precision; K prefixes and all derived views are correlated"},
        "inputs": inputs,
        "saved_old_result": {"path": str(previous_path.relative_to(ROOT)), "sha256": digest(previous_path),
                             "U_source_delete_one_vectors_reused": True},
        "unmarked_inputs": {n: complements[n]["provenance"] for n in NS},
        "increment_receipts": [{"path": str(path.relative_to(ROOT)), "sha256": digest(path)} for path in receipts],
        "code": [{"path": str(path.relative_to(ROOT)), "sha256": digest(path)} for path in source_paths],
        "uncertainty": "Joint first-order covariance with true source/complement batch weights; shared cyclic four-N source group and independent endpoint seed groups; the entire newly marked1M is excluded from each endpoint anchor",
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "scipy": scipy.__version__, "machine": platform.machine()},
        "elapsed_seconds": time.perf_counter() - started,
        "configuration_replays_by_this_script": 0, "server_actions": 0, "test_suites": []}
    OUTPUT.mkdir(parents=True, exist_ok=True)
    with destination.open("x") as handle:
        handle.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"elapsed_seconds": result["elapsed_seconds"], "model_extension_diagnostics": models,
        "source_responses": {n: estimates[f"N{n}.Udot_fugacity"] for n in NS},
        "root_comoving_rank1": {n: estimates[f"N{n}.root_comoving_rank1_fugacity"] for n in NS}}, indent=2))


if __name__ == "__main__":
    main()
