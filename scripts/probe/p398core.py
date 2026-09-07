#!/usr/bin/env python3
"""Independent reimplementation of the P398 exact finite object.

P398 works on the canonical *circular noncrossing connectivity states* of width
``w``: set partitions of the cyclic set {0,...,w-1} with no crossing blocks,
written as restricted-growth strings (RGS).  The number of such states is the
Catalan number C_w.  Every boundary point carries two competing moves:

    join@p   : merge point p with its next cyclic neighbour (p+1) mod w
    detach@p : split point p away from its block into its own singleton

and the declared generator pencil is

    G_eta = (1+eta)*sum(join) + (1-eta)*sum(detach) - exit-rate diagonal
          = G0 + eta*H,      G0 = J + D,  H = J - D.

This module is written from the public problem statement (repository
LightChainr/Matching-One, issues #398/#580/#588/#598) as an independent
implementation, so the later theorem checks do not share code with the
repository scripts.

Conventions kept identical to the repository's declared object:

* states enumerated as canonical RGS in lexicographic order;
* ``wrap`` readout = [ state[0] == state[w-1] ];
* the surviving reflection is the one fixing the cut between w-1 and 0,
  ``R : i -> (w-1) - i``;
* sources: delta on all-singletons, all-in-one-block, wrapped pair, uniform.

Everything below is deterministic and uses integer/fraction arithmetic where a
quantity is claimed to be exact.
"""

from __future__ import annotations

import itertools
from functools import lru_cache
from fractions import Fraction
from typing import Dict, List, Optional, Sequence, Tuple


# --------------------------------------------------------------------------
# 1. state space
# --------------------------------------------------------------------------


def canonical_rgs(labels: Sequence[int]) -> Tuple[int, ...]:
    mapping: Dict[int, int] = {}
    out: List[int] = []
    for label in labels:
        if label not in mapping:
            mapping[label] = len(mapping)
        out.append(mapping[label])
    return tuple(out)


def is_noncrossing(state: Sequence[int]) -> bool:
    """No two blocks interleave: no a<b<c<d with a,c in one block, b,d in another."""
    w = len(state)
    for a in range(w):
        for b in range(a + 1, w):
            if state[a] == state[b]:
                continue
            for c in range(b + 1, w):
                if state[c] != state[a]:
                    continue
                for d in range(c + 1, w):
                    if state[d] == state[b]:
                        return False
    return True


def _all_partitions_rgs(width: int) -> List[Tuple[int, ...]]:
    """All set partitions of {0..w-1} as RGS, in lexicographic order."""
    out: List[Tuple[int, ...]] = []

    def visit(prefix: Tuple[int, ...], maximum: int) -> None:
        if len(prefix) == width:
            out.append(prefix)
            return
        for label in range(maximum + 2):
            visit(prefix + (label,), max(maximum, label))

    visit((0,), 0)
    return out


@lru_cache(maxsize=None)
def states(width: int) -> Tuple[Tuple[int, ...], ...]:
    return tuple(s for s in _all_partitions_rgs(width) if is_noncrossing(s))


def catalan_number(w: int) -> int:
    import math

    return math.comb(2 * w, w) // (w + 1)


def rgs_to_blocks(state: Sequence[int]) -> Tuple[frozenset[int], ...]:
    out = []
    for b in range(max(state) + 1):
        out.append(frozenset(i for i, x in enumerate(state) if x == b))
    return tuple(out)


# --------------------------------------------------------------------------
# 2. moves and generator rows
# --------------------------------------------------------------------------


def join_adjacent(state: Sequence[int], point: int) -> Tuple[int, ...]:
    """Merge ``point`` with its next cyclic neighbour ``(point+1) % w``."""
    w = len(state)
    neighbour = (point + 1) % w
    source = state[neighbour]
    target = state[point]
    labels = [target if x == source else x for x in state]
    return canonical_rgs(labels)


def detach_point(state: Sequence[int], point: int) -> Tuple[int, ...]:
    """Split ``point`` into its own singleton block."""
    labels = list(state)
    labels[point] = max(state) + 1
    return canonical_rgs(labels)


class Generator:
    """Move targets and row-oriented generator for one width.

    ``target[(op, p)][s]`` is the state reached by applying move ``(op, p)`` to
    state index ``s``, or ``None`` when the move does not change the state
    (self transitions are dropped: in a continuous-time chain they cancel
    between the off-diagonal entry and the exit rate).
    """

    def __init__(self, width: int) -> None:
        self.width = width
        self.states = states(width)
        self.size = len(self.states)
        self.rank = {s: i for i, s in enumerate(self.states)}
        self.moves: Tuple[Tuple[str, int], ...] = tuple(
            (op, p) for p in range(width) for op in ("join", "detach")
        )
        self.target: Dict[Tuple[str, int], List[Optional[int]]] = {}
        for p in range(width):
            for op, impl in (
                ("join", join_adjacent),
                ("detach", detach_point),
            ):
                col: List[Optional[int]] = []
                for s in self.states:
                    image = self.rank[impl(s, p)]
                    col.append(None if image == self.rank[s] else image)
                self.target[(op, p)] = col

    def rows(self, rates: Dict[Tuple[str, int], float]) -> List[List[Tuple[int, float]]]:
        """Row-oriented generator; ``row[s]`` = sorted (col, value) pairs."""
        out: List[List[Tuple[int, float]]] = []
        for source in range(self.size):
            entries: Dict[int, float] = {}
            diagonal = 0.0
            for move, rate in rates.items():
                image = self.target[move][source]
                if image is None or rate == 0.0:
                    continue
                entries[image] = entries.get(image, 0.0) + rate
                diagonal -= rate
            entries[source] = entries.get(source, 0.0) + diagonal
            out.append(sorted(entries.items()))
        return out

    def rows_exact(
        self, rates: Dict[Tuple[str, int], Fraction]
    ) -> List[List[Tuple[int, Fraction]]]:
        out: List[List[Tuple[int, Fraction]]] = []
        for source in range(self.size):
            entries: Dict[int, Fraction] = {}
            diagonal = Fraction(0)
            for move, rate in rates.items():
                image = self.target[move][source]
                if image is None or rate == 0:
                    continue
                entries[image] = entries.get(image, Fraction(0)) + rate
                diagonal -= rate
            if diagonal != 0:
                entries[source] = entries.get(source, Fraction(0)) + diagonal
            out.append(sorted(entries.items()))
        return out

    # ---- convenience rate dictionaries -----------------------------------
    def rate_baseline(self) -> Dict[Tuple[str, int], float]:
        return {m: 1.0 for m in self.moves}

    def rate_join(self) -> Dict[Tuple[str, int], float]:
        return {m: (1.0 if m[0] == "join" else 0.0) for m in self.moves}

    def rate_detach(self) -> Dict[Tuple[str, int], float]:
        return {m: (1.0 if m[0] == "detach" else 0.0) for m in self.moves}

    def rate_single_join(self, p: int) -> Dict[Tuple[str, int], float]:
        return {m: (1.0 if m == ("join", p) else 0.0) for m in self.moves}


# --------------------------------------------------------------------------
# 3. readouts and sources
# --------------------------------------------------------------------------


def readout_functions(width: int):
    """The eight declared readouts.  Returns dict name -> function(state)."""
    half = width // 2

    def observable_blocks(state):
        return float(max(state) + 1)

    def observable_singletons(state):
        return float(sum(1 for b in rgs_to_blocks(state) if len(b) == 1))

    def observable_wrap(state):
        return 1.0 if state[0] == state[-1] else 0.0

    def observable_max_block(state):
        return float(max(len(b) for b in rgs_to_blocks(state)))

    def observable_linked_pairs(state):
        return float(
            sum(len(b) * (len(b) - 1) // 2 for b in rgs_to_blocks(state))
        )

    def observable_halves_linked(state):
        for b in rgs_to_blocks(state):
            if any(p < half for p in b) and any(p >= half for p in b):
                return 1.0
        return 0.0

    def observable_covering_depth(state):
        blocks = [sorted(b) for b in rgs_to_blocks(state)]
        best = 0
        for point in range(width):
            depth = sum(
                1
                for b in blocks
                if b[0] < point < b[-1] and point not in b
            )
            best = max(best, depth)
        return float(best)

    def observable_boundary_span(state):
        return float(sum(max(b) - min(b) for b in rgs_to_blocks(state)))

    return {
        "blocks": observable_blocks,
        "singletons": observable_singletons,
        "wrap": observable_wrap,
        "max_block": observable_max_block,
        "linked_pairs": observable_linked_pairs,
        "halves_linked": observable_halves_linked,
        "covering_depth": observable_covering_depth,
        "boundary_span": observable_boundary_span,
    }


PRIMARY = ("blocks", "singletons", "wrap")
HELD_OUT = ("max_block", "linked_pairs", "halves_linked", "covering_depth", "boundary_span")
ALL_READOUTS = PRIMARY + HELD_OUT


def source_vectors(g: Generator) -> Dict[str, List[float]]:
    """The four frozen initial laws (delta distributions + uniform)."""
    r = {s: i for i, s in enumerate(g.states)}
    w = g.width
    v: Dict[str, List[float]] = {}
    names = (
        ("delta_all_singletons", tuple(range(w))),
        ("delta_single_block", tuple(0 for _ in range(w))),
        ("delta_wrapped_pair", tuple([0] + list(range(1, w - 1)) + [0])),
    )
    for name, state in names:
        vec = [0.0] * g.size
        vec[r[state]] = 1.0
        v[name] = vec
    v["uniform"] = [1.0 / g.size] * g.size
    return v


# --------------------------------------------------------------------------
# 4. reflection and parity
# --------------------------------------------------------------------------


def reflection_perm(width: int, k: int) -> List[int]:
    """Reflection i -> (k - i) mod w acting on point labels."""
    return [(k - i) % width for i in range(width)]


def apply_permutation(state: Sequence[int], perm: Sequence[int]) -> Tuple[int, ...]:
    """State obtained by moving point i to position perm[i]."""
    width = len(state)
    inverse = [0] * width
    for source, destination in enumerate(perm):
        inverse[destination] = source
    return canonical_rgs(state[inverse[j]] for j in range(width))


def state_permutation(g: Generator, perm: Sequence[int]) -> List[int]:
    """Induced permutation of the *state index list*."""
    r = {s: i for i, s in enumerate(g.states)}
    return [r[apply_permutation(s, perm)] for s in g.states]


def orbit_labels(g: Generator, pi: Sequence[int]) -> List[int]:
    labels = [-1] * g.size
    next_label = 0
    for i in range(g.size):
        if labels[i] < 0:
            labels[i] = labels[pi[i]] = next_label
            next_label += 1
    return labels


def matvec(rows, vector):
    return [sum(value * vector[col] for col, value in row) for row in rows]


def rows_operator(rows, vec):
    """rows act on column functions from the left: (G f)_i = sum_j G_ij f_j."""
    return matvec(rows, vec)
