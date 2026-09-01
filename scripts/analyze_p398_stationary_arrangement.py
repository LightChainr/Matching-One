#!/usr/bin/env python3
"""Locate stationary response in block sizes, dual sizes and cycle arrangement."""
from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import io
import json
from pathlib import Path
import platform
import subprocess
import sys
import time
import types

import numpy as np
import scipy
from scipy import linalg

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis/p398_stationary_arrangement_contract.json"
OUTPUT = ROOT / "results/p398-stationary-arrangement/latest.json"


def encode(value):
    a = np.asarray(value)
    return np.stack((a.real, a.imag), axis=-1).tolist()


def decode(value):
    a = np.asarray(value)
    return a[..., 0] + 1j * a[..., 1]


def maximum(value):
    return float(np.max(np.abs(value)))


def main():
    started = time.perf_counter()
    if OUTPUT.exists():
        raise ValueError("Preserve the saved result; reproduce in another output copy")
    contract = json.loads(CONTRACT.read_text())
    inputs = []

    def read(commit, path):
        data = subprocess.check_output(["git", "show", f"{commit}:{path}"], cwd=ROOT)
        inputs.append({"commit": commit, "path": path,
                       "sha256": hashlib.sha256(data).hexdigest()})
        return data

    def source(name):
        return read(contract["source_commit"], contract["source_directory"] + "/" + name)

    # Source modules only supply unchanged state ordering and generator geometry.
    for name in ("frozen_model", "frozen_rate_model"):
        code = source(name + ".py")
        module = types.ModuleType(name)
        module.__file__ = f"{contract['source_commit']}/{name}.py"
        sys.modules[name] = module
        exec(compile(code, module.__file__, "exec"), module.__dict__)
    codec, model = sys.modules["frozen_model"], sys.modules["frozen_rate_model"]
    saved = json.loads(read(contract["score_commit"], contract["score_path"]))
    old = json.loads(source("results/latest.json"))
    archived = np.load(io.BytesIO(source("previous_character_i_zero.npz")))
    states, _, q, rays, _, _, complement = model.construct(contract["width"])
    pi = np.asarray(saved["stationary_measure"]["pi"])
    score = np.asarray(saved["stationary_measure"]["true_score"])
    if len(states) != len(pi) or maximum(rays - archived["source"]) > 1e-12:
        raise ValueError("source state/basis ordering does not match the saved score")
    state_index = {s: i for i, s in enumerate(states)}
    profiles = [tuple(sorted(Counter(s).values(), reverse=True)) for s in states]
    dual_profiles = [profiles[i] for i in complement]
    joint_profiles = list(zip(profiles, dual_profiles))

    def orbit_key(state):
        orbit = []
        for seq in (state, tuple(reversed(state))):
            for step in range(len(state)):
                orbit.append(codec.canonical_rgs(seq[step:] + seq[:step]))
        return min(orbit)

    coordinates = [
        ("block_count", [len(p) for p in profiles]),
        ("primal_sizes", profiles),
        ("primal_dual_sizes", joint_profiles),
        ("dihedral_arrangement", [orbit_key(s) for s in states]),
    ]
    variance = float(pi @ (score - pi @ score)**2)
    h = archived["mass"]
    inverse = linalg.inv(h)
    c0 = decode(old["integrals"]["C0"])
    integrated_u = decode(old["integrals"]["integrated_U"])
    dynamic = decode(old["integrals"]["integrated_Uprime_generator"])

    def response(g):
        metric_prime = (q.conj().T @ q.multiply((pi*g)[:, None])).toarray()
        c0prime = (rays.conj().T @ metric_prime @ rays).conj()
        integrated_cprime = (rays.conj().T @ metric_prime @ inverse @ rays).conj()
        stationary = linalg.solve(c0, integrated_cprime - c0prime @ integrated_u)
        return c0prime, stationary

    def state_record(i):
        return {"index": int(i), "partition": list(states[i]), "pi": float(pi[i]),
                "score": float(score[i]), "primal_sizes": list(profiles[i]),
                "dual_sizes": list(dual_profiles[i])}

    stages, components = [], []
    previous = np.full_like(score, float(pi @ score))
    previous_key = None
    joint_groups = None
    for name, keys in coordinates:
        groups = defaultdict(list)
        for i, key in enumerate(keys):
            groups[key].append(i)
        projected = np.zeros_like(score)
        group_rows, witness = [], None
        for key, members in sorted(groups.items()):
            ix = np.asarray(members)
            if previous_key is not None and len({previous_key[i] for i in members}) != 1:
                raise ValueError("coordinate hierarchy is not nested")
            mass = float(pi[ix].sum())
            projected[ix] = float(pi[ix] @ score[ix] / mass)
            lo, hi = ix[np.argmin(score[ix])], ix[np.argmax(score[ix])]
            span = float(score[hi] - score[lo])
            group_rows.append({"coordinate": key, "states": len(members), "pi_mass": mass,
                               "mean_score": float(projected[lo]), "score_span": span})
            if witness is None or span > witness["score_span"]:
                witness = {"coordinate": key, "score_span": span,
                           "states": [state_record(lo), state_record(hi)]}
        residual = float(pi @ (score - projected)**2)
        increment = projected - previous
        contribution = float(pi @ increment**2)
        cprime, stationary = response(projected)
        dc, dr = response(increment)
        stages.append({"id": name, "group_count": len(groups),
                       "residual_variance": residual, "residual_fraction": residual/variance,
                       "added_explained_variance": contribution,
                       "added_explained_fraction": contribution/variance,
                       "Kreweras_odd_max_residual": maximum(projected[complement] + projected),
                       "static_C0prime": encode(cprime),
                       "stationary_integrated_Uprime": encode(stationary),
                       "total_integrated_Uprime_with_saved_generator": encode(stationary + dynamic),
                       "projected_score": projected.tolist(), "groups": group_rows,
                       "same_coordinate_witness": witness})
        components.append({"id": name, "score_variance": contribution,
                           "static_C0prime": encode(dc),
                           "stationary_integrated_Uprime": encode(dr)})
        previous, previous_key = projected, keys
        if name == "primal_dual_sizes":
            joint_groups = groups

    remainder_c, remainder_r = response(score - previous)
    components.append({"id": "within_dihedral_remainder",
                       "score_variance": float(pi @ (score - previous)**2),
                       "static_C0prime": encode(remainder_c),
                       "stationary_integrated_Uprime": encode(remainder_r)})
    true_c, true_r = response(score)
    mean_c, mean_r = response(np.full_like(score, float(pi @ score)))
    closure_r = sum((decode(c["stationary_integrated_Uprime"]) for c in components), mean_r)

    # Exact rate counts test the physical size-profile quotient, not software.
    group_keys = sorted(joint_groups)
    group_index = {key: i for i, key in enumerate(group_keys)}
    counts = np.zeros((len(states), 2, len(group_keys)), dtype=np.int64)
    for i, state in enumerate(states):
        for channel, action in enumerate((model.join_adjacent, model.detach_state)):
            for site in range(contract["width"]):
                destination = state_index[action(state, site)]
                counts[i, channel, group_index[joint_profiles[destination]]] += 1
    lumpability = {"profile_classes": len(group_keys), "attempts_per_channel": contract["width"]}
    for label, vectors in (("baseline_J_plus_D", counts[:, 0] + counts[:, 1]),
                           ("tangent_J_minus_D", counts[:, 0] - counts[:, 1])):
        failed, witness = 0, None
        for key in group_keys:
            members = joint_groups[key]
            reference = members[0]
            differing = next((i for i in members if np.any(vectors[i] != vectors[reference])), None)
            if differing is None:
                continue
            failed += 1
            if witness is None:
                target = int(np.flatnonzero(vectors[differing] != vectors[reference])[0])
                pair = [reference, differing]
                witness = {"source_profile": key, "target_profile": group_keys[target],
                           "states": [state_record(i) for i in pair],
                           "join_attempts_to_target": [int(counts[i, 0, target]) for i in pair],
                           "detach_attempts_to_target": [int(counts[i, 1, target]) for i in pair],
                           "compared_rates": [int(vectors[i, target]) for i in pair],
                           "self_class": target == group_index[key],
                           "diagonal_note": "Each source has the same -8 per-channel diagonal subtraction, so count differences equal generator differences even in the source class."}
        lumpability[label] = {"strongly_lumpable": failed == 0,
                             "nonconstant_rate_profile_classes": failed, "witness": witness}

    result = {
        "schema": "matching-one.p398-stationary-arrangement.v1",
        "status": "computed_conditional_stationary_score_and_exact_profile_transition_counts",
        "contract": contract,
        "execution_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "contract_sha256": hashlib.sha256(CONTRACT.read_bytes()).hexdigest(),
        "inputs": inputs, "state_count": len(states), "character_dimension": len(h),
        "score_mean": float(pi @ score), "score_variance": variance,
        "state_records": [state_record(i) for i in range(len(states))],
        "stages": stages, "components": components,
        "true_static_C0prime": encode(true_c), "true_stationary_integrated_Uprime": encode(true_r),
        "reused_true_generator_term": encode(dynamic),
        "true_total_integrated_Uprime": old["integrals"]["integrated_Uprime"],
        "score_mean_static_roundoff": encode(mean_c),
        "score_mean_integral_roundoff": encode(mean_r),
        "additive_integral_max_residual": maximum(closure_r - true_r),
        "saved_integral_max_residual": maximum(true_r + dynamic - decode(old["integrals"]["integrated_Uprime"])),
        "explained_variance_addback_residual": float(sum(c["score_variance"] for c in components) - variance),
        "profile_lumpability": lumpability,
        "environment": {"python": platform.python_version(), "platform": platform.platform(),
                        "numpy": np.__version__, "scipy": scipy.__version__},
        "elapsed_seconds": time.perf_counter() - started,
        "stationary_solves": 0, "new_eta_points": 0, "Frechet_recomputations": 0,
        "new_MC_samples": 0, "server_operations": 0, "test_suites": 0,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"elapsed_seconds": result["elapsed_seconds"],
                      "stages": [{k: r[k] for k in ("id", "group_count", "residual_fraction", "added_explained_fraction", "Kreweras_odd_max_residual", "same_coordinate_witness", "total_integrated_Uprime_with_saved_generator")} for r in stages],
                      "lumpability": lumpability,
                      "integral_closure": result["additive_integral_max_residual"]}, indent=2))


if __name__ == "__main__":
    main()
