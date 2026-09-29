#!/usr/bin/env python3
"""Exact L3 occupied-set geometry and first-birth history transport.

512 subsets per lattice, dynamic integer prefix counts; no 9! permutation
rerun and no Monte Carlo. New scientific readout: completion-pivotal count
and primitive rank-one homology direction.
"""
from collections import Counter, defaultdict
from fractions import Fraction
import gzip
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path
import platform
import time


HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent / "birth-gap-20260929/data"
L, N = 3, 9
FACT = [math.factorial(k) for k in range(N+1)]


def primitive(x, y):
    g = math.gcd(x, y)
    x, y = x//g, y//g
    if x < 0 or (x == 0 and y < 0):
        x, y = -x, -y
    return x, y


def topology(mask, neighbors):
    """Lifted graph traversal; ambient image span, not max component rank."""
    potentials = {}
    basis = None
    rank = 0
    for root in range(N):
        if not mask >> root & 1 or root in potentials:
            continue
        potentials[root] = (0, 0)
        stack = [root]
        while stack:
            u = stack.pop()
            ux, uy = potentials[u]
            for v, dx, dy in neighbors[u]:
                if not mask >> v & 1:
                    continue
                expected = ux+dx, uy+dy
                if v not in potentials:
                    potentials[v] = expected
                    stack.append(v)
                else:
                    vx, vy = potentials[v]
                    x, y = expected[0]-vx, expected[1]-vy
                    assert x % L == 0 and y % L == 0
                    x, y = x//L, y//L
                    if x or y:
                        if basis is None:
                            basis = primitive(x, y)
                            rank = 1
                        elif basis[0]*y != basis[1]*x:
                            rank = 2
    return rank, basis if rank == 1 else None


def bits(mask):
    return [i for i in range(N) if mask >> i & 1]


def analyze(lattice):
    steps = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    if lattice == "triangular":
        steps += [(1, 1), (-1, -1)]
    neighbors = [[(((v % L + dx) % L)+L*((v//L+dy) % L), dx, dy)
                  for dx, dy in steps] for v in range(N)]
    topo = [topology(mask, neighbors) for mask in range(1 << N)]
    ways = [Counter() for _ in topo]
    rows, hist = [], Counter()
    rank_counts = Counter()
    directions = set()
    # Proper subsets have smaller integer masks: this is a topological DP order.
    for mask, (rank, direction) in enumerate(topo):
        k = mask.bit_count()
        rank_counts[k, rank] += 1
        if rank == 1:
            directions.add(direction)
            for v in bits(mask):
                old = mask ^ (1 << v)
                previous_rank, previous_direction = topo[old]
                if previous_rank == 0:
                    ways[mask][k] += FACT[k-1]
                else:
                    assert previous_rank == 1 and previous_direction == direction
                    ways[mask].update(ways[old])
            assert sum(ways[mask].values()) == FACT[k]
            pivotal = [v for v in range(N) if not mask >> v & 1
                       and topo[mask | (1 << v)][0] == 2]
            nu = len(pivotal)
            for j1, count in ways[mask].items():
                hist[j1, k+1, direction] += count*nu*FACT[N-k-1]
            rows.append({"mask": mask, "occupied": bits(mask), "k": k,
                         "direction": list(direction), "completion_sites": pivotal,
                         "nu2": nu, "prefix_birth_counts": dict(sorted(ways[mask].items()))})
        elif rank == 0 and k < N:
            direct = sum(topo[mask | (1 << v)][0] == 2
                         for v in range(N) if not mask >> v & 1)
            hist[k+1, k+1, None] += FACT[k]*direct*FACT[N-k-1]
    hist = Counter({key: value for key, value in hist.items() if value})
    assert sum(hist.values()) == FACT[N]
    # A single narrow source-convention control, not a repeat of acquisition.
    filename = f"exact-control-{lattice}-L3-b00.json.gz"
    raw = (ARCHIVE / filename).read_bytes()
    archived = json.loads(gzip.decompress(raw))
    collapsed = Counter()
    for (j1, j2, direction), count in hist.items():
        collapsed[j1, j2] += count
    assert collapsed == Counter({(j1, j2): count for j1, j2, count in archived["histogram"]})
    assert archived["samples"] == FACT[N]

    profiles = []
    for k in range(N):
        at_k = [r for r in rows if r["k"] == k]
        if not at_k:
            continue
        counts = Counter(r["nu2"] for r in at_k)
        profiles.append({"k": k, "rank_one_configurations": len(at_k),
                         "nu2_state_counts": dict(sorted(counts.items()))})

    direction_kernels = []
    for direction in sorted(directions):
        subhist = {(j1, j2): count for (j1, j2, d), count in hist.items() if d == direction}
        h = lambda s, t: sum(count for (j1, j2), count in subhist.items() if j1 <= s and j2 > t)
        violations = []
        for a, b, c in combinations(range(N+1), 3):
            A, B, C, D = h(a, b), h(a, c), h(b, c), h(b, b)
            det = B*D-A*C
            if det:
                violations.append({"times": [a, b, c],
                                   "conditional_covariance": str(Fraction(det, D*D)),
                                   "survival_difference": str(Fraction(B, A)-Fraction(C-B, D-A))})
        direction_kernels.append({
            "direction": list(direction), "non_direct_permutation_count": sum(subhist.values()),
            "paired_histogram": [[j1, j2, count] for (j1, j2), count in sorted(subhist.items())],
            "count_time_violations": violations,
        })

    # Exact decomposition of the previously used one-step witness at (4,5,6).
    risk_rows = [r for r in rows if r["k"] == 5]
    by_direction_and_nu = defaultdict(lambda: [0, 0])
    for r in risk_rows:
        early = sum(count for j1, count in r["prefix_birth_counts"].items() if j1 <= 4)
        late = r["prefix_birth_counts"].get(5, 0)
        by_direction_and_nu[tuple(r["direction"]), r["nu2"]][0] += early
        by_direction_and_nu[tuple(r["direction"]), r["nu2"]][1] += late
    early_total = sum(x[0] for x in by_direction_and_nu.values())
    late_total = sum(x[1] for x in by_direction_and_nu.values())
    early_mean = Fraction(sum(nu*c[0] for (d, nu), c in by_direction_and_nu.items()), early_total)
    late_mean = Fraction(sum(nu*c[1] for (d, nu), c in by_direction_and_nu.items()), late_total)
    contrast = -(early_mean-late_mean)/(N-5)

    # Strong-lumping diagnostic on rank-one rows only. It concerns arbitrary
    # preparations, not automatically the weak Markov property from empty.
    coarse = defaultdict(list)
    row_lookup = {r["mask"]: r for r in rows}
    for r in rows:
        signature = Counter()
        mask = r["mask"]
        for v in range(N):
            if mask >> v & 1:
                continue
            new = mask | (1 << v)
            if topo[new][0] == 2:
                signature["absorbed_rank2"] += 1
            else:
                next_row = row_lookup[new]
                signature[f"rank1_nu{next_row['nu2']}"] += 1
        coarse[r["k"], tuple(r["direction"]), r["nu2"]].append((mask, tuple(sorted(signature.items()))))
    failures = []
    for key, members in sorted(coarse.items()):
        distinct = {}
        for mask, signature in members:
            distinct.setdefault(signature, mask)
        if len(distinct) > 1:
            failures.append({"k": key[0], "direction": list(key[1]), "nu2": key[2],
                             "witnesses": [{"mask": mask, "successor_counts": dict(sig)}
                                           for sig, mask in list(distinct.items())[:2]]})
    # Drop the direction and retain only one completion phase at each k.
    # These rows establish strong closure within rank one, not within rank zero.
    by_k_nu = defaultdict(list)
    for (k, direction, nu), members in coarse.items():
        by_k_nu[k, nu].extend(members)
    transition_classes = []
    for (k, nu), members in sorted(by_k_nu.items()):
        signatures = {sig for mask, sig in members}
        assert len(signatures) == 1
        signature = dict(next(iter(signatures)))
        transition_classes.append({"k": k, "nu2": nu, "states": len(members),
                                   "successor_counts": signature,
                                   "denominator": N-k})

    # An exact count-clock weak Markov realization from the empty start.
    # States are zero, rank-one slow, rank-one fast, two; k is external time.
    min_nu = {k: min(nu for kk, nu in by_k_nu if kk == k) for k, nu in by_k_nu}
    def phase(mask):
        rank = topo[mask][0]
        if rank == 0:
            return 0
        if rank == 2:
            return 3
        r = row_lookup[mask]
        return 1+int(r["nu2"] != min_nu[r["k"]])
    matrices = []
    for k in range(N):
        classes = defaultdict(list)
        for mask in range(1 << N):
            if mask.bit_count() == k:
                classes[phase(mask)].append(mask)
        matrix = [[Fraction(int(i == j)) for j in range(4)] for i in range(4)]
        reachable = []
        for state, masks in classes.items():
            reachable.append(state)
            if state == 3:
                continue
            counts = Counter()
            signatures = set()
            for mask in masks:
                local = Counter(phase(mask | (1 << v)) for v in range(N) if not mask >> v & 1)
                counts.update(local)
                signatures.add(tuple(sorted(local.items())))
            if state in (1, 2):
                assert len(signatures) == 1
            matrix[state] = [Fraction(counts[j], len(masks)*(N-k)) for j in range(4)]
        assert all(sum(row) == 1 for row in matrix)
        matrices.append({"k": k, "reachable_states": sorted(reachable),
                         "matrix": [[str(x) for x in row] for row in matrix]})

    # Reproduce the full paired-birth law using this small chain, without paths.
    law = {(0, 0): Fraction(1)}
    chain_hist = Counter()
    for entry in matrices:
        k = entry["k"]
        matrix = [[Fraction(x) for x in row] for row in entry["matrix"]]
        new_law = defaultdict(Fraction)
        for (state, j1), mass in law.items():
            for next_state, probability in enumerate(matrix[state]):
                if not probability:
                    continue
                next_j1 = k+1 if state == 0 and next_state != 0 else j1
                if next_state == 3:
                    chain_hist[next_j1, k+1] += mass*probability
                else:
                    new_law[next_state, next_j1] += mass*probability
        law = dict(new_law)
    assert not law
    assert chain_hist == Counter({key: Fraction(value, FACT[N]) for key, value in collapsed.items()})
    return {
        "lattice": lattice, "L": L, "N": N, "subsets": 1 << N,
        "archived_histogram_matches": True,
        "archive_file": "analysis/birth-gap-20260929/data/"+filename,
        "archive_sha256": hashlib.sha256(raw).hexdigest(),
        "primitive_unoriented_directions": [list(d) for d in sorted(directions)],
        "rank_counts": [[k, rank, count] for (k, rank), count in sorted(rank_counts.items())],
        "completion_profiles": profiles,
        "direction_kernels": direction_kernels,
        "risk_k5_entry_cut4": {
            "grain": "ordered prefixes ending at rank-one occupied set; multiply by 4! for full-permutation counts",
            "early_prefixes": early_total, "late_prefixes": late_total,
            "mean_nu2_early": str(early_mean), "mean_nu2_late": str(late_mean),
            "one_step_survival_difference": str(contrast),
            "composition": [{"direction": list(d), "nu2": nu, "early_prefixes": counts[0], "late_prefixes": counts[1]}
                            for (d, nu), counts in sorted(by_direction_and_nu.items())],
        },
        "rank_one_k_direction_nu2_strong_lumping_failures": failures,
        "rank_one_k_nu2_transition_classes": transition_classes,
        "count_clock_phase_realization": {
            "states": ["rank0", "rank1_slow", "rank1_fast", "rank2"],
            "phase": "at fixed k, slow is minimum nu2 among rank-one states; fast is the other value if present",
            "external_time": "occupied count k, not hidden iid-label count",
            "claim": "weak Markov realization from the empty uniform-permutation start; rank-one rows strongly close across all configurations",
            "off_risk_rows": "identity placeholders; not a statement about artificial unreachable preparations",
            "paired_birth_histogram_reproduced_exactly": True,
            "transition_matrices": matrices,
        },
        "rank_one_configuration_rows": rows,
    }


def main():
    started = time.perf_counter()
    cases = [analyze(lattice) for lattice in ("square", "triangular")]
    report = {
        "schema": "matching-one.birth-completion-geometry.v1",
        "as_of": "2026-09-29",
        "question": "Does rank-one direction remove birth-history memory, and what immediate completion geometry carries it?",
        "method": "lifted traversal of 512 subsets per lattice; exact prefix dynamic programming and outgoing-site counts",
        "direction": "primitive unoriented rational homology line; fixed until rank two by image inclusion",
        "nu2": "number of vacant vertices whose addition changes ambient rank from one to two",
        "one_step_exit": "nu2/(N-k) conditional on the full occupied configuration at count k",
        "cases": cases,
        "no_monte_carlo_or_permutation_enumeration": True,
        "limitations": ["L3 only; not large-size or continuum sufficiency.",
                        "One-step hazard sufficiency is not recursive state sufficiency.",
                        "Strong-lumping failure is not by itself failure of the weak Markov law from the empty state.",
                        "A homology direction is not the full connectivity or geometric state."],
        "python": platform.python_version(), "platform": platform.platform(),
        "wall_seconds": time.perf_counter()-started,
    }
    (HERE / "result.json").write_text(json.dumps(report, indent=2)+"\n")
    for c in cases:
        print(json.dumps({k: c[k] for k in ("lattice", "primitive_unoriented_directions", "completion_profiles",
                                           "risk_k5_entry_cut4", "rank_one_k_direction_nu2_strong_lumping_failures")}, indent=2))
        print("direction violation counts",[len(d["count_time_violations"]) for d in c["direction_kernels"]])
    print("seconds", report["wall_seconds"])


if __name__ == "__main__":
    main()
