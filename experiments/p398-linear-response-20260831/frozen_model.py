"""Unchanged function bodies extracted from frozen Matching-One sources."""
from __future__ import annotations
from functools import lru_cache
from typing import Iterable, Sequence
import numpy as np


# 1f19fc1a2d9fc59dce650e95268c716762725985:scripts/noncrossing_connectivity_codec.py:19
def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


# 1f19fc1a2d9fc59dce650e95268c716762725985:scripts/noncrossing_connectivity_codec.py:24
def canonical_rgs(labels: Iterable[int]) -> tuple[int, ...]:
    """Canonicalize arbitrary integer block labels by first occurrence."""

    mapping: dict[int, int] = {}
    output: list[int] = []
    for label in labels:
        _require(isinstance(label, int), "block labels must be integers")
        if label not in mapping:
            mapping[label] = len(mapping)
        output.append(mapping[label])
    return tuple(output)


# 1f19fc1a2d9fc59dce650e95268c716762725985:scripts/noncrossing_connectivity_codec.py:37
def validate_rgs(state: Sequence[int]) -> tuple[int, ...]:
    state = tuple(state)
    _require(bool(state), "connectivity state must not be empty")
    _require(all(isinstance(label, int) for label in state), "block labels must be integers")
    _require(state[0] == 0, "restricted-growth string must start at zero")
    maximum = 0
    for index, label in enumerate(state[1:], start=1):
        _require(0 <= label <= maximum + 1, f"invalid restricted-growth label at index {index}")
        maximum = max(maximum, label)
    _require(canonical_rgs(state) == state, "restricted-growth string is not canonical")
    return state


# 1f19fc1a2d9fc59dce650e95268c716762725985:scripts/noncrossing_connectivity_codec.py:50
def is_noncrossing_rgs(state: Sequence[int]) -> bool:
    """Return false exactly when two blocks contain an alternating a<b<c<d."""

    state = validate_rgs(state)
    width = len(state)
    for a in range(width):
        for b in range(a + 1, width):
            if state[a] == state[b]:
                continue
            for c in range(b + 1, width):
                if state[c] != state[a]:
                    continue
                for d in range(c + 1, width):
                    if state[d] == state[b]:
                        return False
    return True


# 1f19fc1a2d9fc59dce650e95268c716762725985:scripts/noncrossing_connectivity_codec.py:68
def generate_rgs_partitions(width: int) -> tuple[tuple[int, ...], ...]:
    """Generate set partitions directly as restricted-growth strings."""

    _require(width >= 1, "width must be positive")
    output: list[tuple[int, ...]] = []

    def visit(prefix: tuple[int, ...], maximum: int) -> None:
        if len(prefix) == width:
            output.append(prefix)
            return
        for label in range(maximum + 2):
            visit(prefix + (label,), max(maximum, label))

    visit((0,), 0)
    return tuple(output)


# 1f19fc1a2d9fc59dce650e95268c716762725985:scripts/noncrossing_connectivity_codec.py:132
@lru_cache(maxsize=None)
def noncrossing_states(width: int) -> tuple[tuple[int, ...], ...]:
    return tuple(state for state in generate_rgs_partitions(width) if is_noncrossing_rgs(state))


# 1f19fc1a2d9fc59dce650e95268c716762725985:scripts/p321_homology_trace_certificate.py:111
def join_adjacent(state: State, site: int) -> State:
    width = len(state)
    other = (site + 1) % width
    left_label = state[site]
    right_label = state[other]
    if left_label == right_label:
        return state
    return canonical_rgs(
        left_label if label == right_label else label for label in state
    )


# 1f19fc1a2d9fc59dce650e95268c716762725985:scripts/p321_homology_trace_certificate.py:103
def rotate_state(state: State, steps: int = 1) -> State:
    width = len(state)
    rotated = [0] * width
    for site, label in enumerate(state):
        rotated[(site + steps) % width] = label
    return canonical_rgs(rotated)


# 1f19fc1a2d9fc59dce650e95268c716762725985:scripts/p333_generic_q_detach_intertwiner.py:47
def detach_state(state: Sequence[int], site: int) -> tuple[int, ...]:
    """Split ``site`` into a singleton; an existing singleton is unchanged."""

    state = tuple(state)
    label = state[site]
    if state.count(label) == 1:
        return state
    fresh = max(state) + 1
    output = list(state)
    output[site] = fresh
    return canonical_rgs(output)


# 1f19fc1a2d9fc59dce650e95268c716762725985:scripts/p398_width8_source_spectrum.py:24
def kreweras(state):
    """Noncrossing complement: cycles of p^-1 c, with c=(0 1 ... w-1)."""
    n = len(state)
    permutation = list(range(n))
    for label in set(state):
        block = [j for j in range(n) if state[j] == label]
        for current, following in zip(block, block[1:]+block[:1]):
            permutation[current] = following
    inverse = [permutation.index(j) for j in range(n)]
    complement = [inverse[(j+1) % n] for j in range(n)]
    labels = [-1]*n
    for j in range(n):
        if labels[j] < 0:
            current = j
            while labels[current] < 0:
                labels[current] = j
                current = complement[current]
    return canonical_rgs(labels)


# 1f19fc1a2d9fc59dce650e95268c716762725985:scripts/p398_width8_memory_motifs.py:65
def next_pair_motifs(states):
    weight=(1j)**np.arange(8)
    output=[]
    for state in states:
        size=[state.count(state[j]) for j in range(8)]
        t3=sum(weight[j]*(size[j]==3) for j in range(8))
        s11=sum((weight[j]+weight[(j+1)%8])*(size[j]==size[(j+1)%8]==1) for j in range(8))
        boundary_two=0j
        for label in set(state):
            members=[j for j in range(8) if state[j]==label]
            if len(members)==2:
                charge=sum(weight[j] for j in members)
                boundary=sum((state[j]==label)!=(state[(j+1)%8]==label) for j in range(8))
                boundary_two+=charge*boundary
        output.append((t3,s11,boundary_two))
    return np.array(output)


# 1f19fc1a2d9fc59dce650e95268c716762725985:scripts/p398_width8_geometric_compression.py:29
def features(states,f,t2):
    pair=next_pair_motifs(states)
    output={"A":f[:,0],"T2":t2,"T3":pair[:,0],"S11":pair[:,1],"B2":pair[:,2]}
    extra=[]
    weight=(1j)**np.arange(8)
    for state in states:
        blocks={label:[j for j in range(8) if state[j]==label] for label in set(state)}
        charge={label:sum(weight[j] for j in sites) for label,sites in blocks.items()}
        sizes={label:len(sites) for label,sites in blocks.items()}
        triplet=[sum(weight[j]*(sizes[state[j]]==3)*
            (sum(state[(j+d)%8]==state[j] for d in (-1,1))==r) for j in range(8)) for r in range(3)]
        t4=sum(charge[label] for label in blocks if sizes[label]==4)
        q3=b3=0j
        for j in range(8):
            a,b=state[j],state[(j+1)%8]
            if a!=b:
                q3+=(charge[a]+charge[b])*(sizes[a]+sizes[b]==3)
                b3+=charge[a]*(sizes[a]==3)+charge[b]*(sizes[b]==3)
        extra.append(triplet+[t4,q3,b3])
    for name,column in zip(("T3_r0","T3_r1","T3_r2","T4","Q3","B3"),np.array(extra).T):
        output[name]=column
    assert np.array_equal(output["T3"],output["T3_r0"]+output["T3_r1"]+output["T3_r2"])
    return output
