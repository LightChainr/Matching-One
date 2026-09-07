#!/usr/bin/env python3
"""Exact gate for the Q-tangent split of issue #581, on tiny square-bond tori.

#581 asks where, in scale and in Boolean degree, a lattice observable's response
to the exact-critical Q score lives: on the duality-**even** Betti tangent, or on
the **ambient homology** source.  Before any Monte Carlo it asks for a
deterministic gate, and this module is that gate.  It adds no evidence; it
establishes that the objects are what they are claimed to be.

**The split.**  On an ``L x L`` square-bond torus at ``p = 1/2``, with ``A`` the
open-edge set, ``A*`` its geometric dual transport, ``k`` the component count and
``b1`` the cycle rank,

    T(A)  = k(A) + |A|/2 = V/2 + B(A)/2,     B = b0 + b1,
    X(A)  = r(A) - 1,     r(A) = rank im[H1(A) -> H1(T^2)],
    T(A) - T(A*) = X(A).

Writing ``B_even = [B(A) + B(A*)]/2``, the score splits as
``T = V/2 + B_even/2 + X/2`` and therefore, for every observable ``O``,

    Cov(O, T) = (1/2) Cov(O, B_even) + (1/2) Cov(O, X),

one duality-even Betti piece and one ambient-homology piece.  Both are verified
here configuration by configuration, not asserted.

**What is checked, and what each check would catch.**

1. ``T - T* = X`` on every configuration.  The other two identities of #581's
   gate, ``T = V/2 + B/2`` and ``B - B* = 2X``, follow algebraically from it and
   are asserted rather than presented as independent evidence.
2. ``r(A)``, the rank of the image in ``H1(T^2)``, computed as the span of every
   component's winding lattice and compared with the largest single component's
   rank.  Two disjoint non-contractible cycles on a torus must be parallel, so
   the two agree -- but that is an argument, and this measures it.
3. The covariance split, exactly, for several observables.
4. A nested spatial filtration with ``Gamma_j = E[D_j D_j^T]``, checked for
   positive semi-definiteness and for telescoping to the full covariance.
5. The Boolean cross-spectrum by Walsh degree, which must sum to the same
   covariance.  This is the ``rho = 1`` end of #227's noise semigroup, where the
   noise operator is the identity.

Everything is exact: at ``p = 1/2`` every configuration is equiprobable, so every
expectation is a dyadic rational and is carried as one.

**Claim boundary.**  Nothing here says a lattice response loads on one piece
rather than the other.  ``B_even`` is a duality-even measure tangent and not
automatically a local energy field; ``X`` is topological but a covariance with it
does not by itself demonstrate a defect theory.  No H4 readout appears, because
the square-*site* matching-odd observable this repository studies has no
canonical lift to a square-bond torus of side 2 or 3, and #581 says explicitly
not to invent one.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

try:  # pragma: no cover - import shape depends on how the script is invoked
    from scripts.square_bond_kappa3 import square_bond_pairs
    from scripts.square_bond_duality_exact import geometric_dual_mask
    from scripts.torus_homology import HomologyUnionFind, _extend_basis
except ModuleNotFoundError:  # pragma: no cover
    from square_bond_kappa3 import square_bond_pairs
    from square_bond_duality_exact import geometric_dual_mask
    from torus_homology import HomologyUnionFind, _extend_basis

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "qtangent-scale-decomposition" / "latest.json"
SCHEMA = "matching-one.qtangent-scale-decomposition.v1"
ISSUE = 581

#: Observables the split is verified against.  Every one is a channel this
#: repository already defines; none is invented here.
OBSERVABLES = ("wrap_either", "wrap_cross", "wrap_direction_0", "open_edges",
               "components", "cycle_rank")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def configuration_statistics(length: int, mask: int, pairs: Sequence[Any]) -> dict[str, int]:
    """Component count, open-edge count, cycle rank and ambient homology rank.

    ``span_rank`` is the rank of the image of ``H1(A)`` in ``H1(T^2)`` -- the
    span of every component's winding lattice.  ``max_component_rank`` is the
    largest single component's, kept so the gate can compare them.
    """
    vertices = length * length
    union_find = HomologyUnionFind(vertices, (length, length))
    edges = 0
    for index, pair in enumerate(pairs):
        if (mask >> index) & 1:
            union_find.add_edge(*pair.primal)
            edges += 1
    roots: set[int] = set()
    basis: list[Any] = []
    max_rank = 0
    for vertex in range(vertices):
        root, _, _ = union_find.find(vertex)
        if root in roots:
            continue
        roots.add(root)
        max_rank = max(max_rank, union_find.component(vertex).rank)
        for winding in union_find.basis[root]:
            _extend_basis(basis, winding)
    components = len(roots)
    return {
        "components": components,
        "open_edges": edges,
        "cycle_rank": edges - vertices + components,
        "span_rank": len(basis),
        "max_component_rank": max_rank,
    }


def wrapping_flags(length: int, mask: int, pairs: Sequence[Any]) -> dict[str, int]:
    """The repository's wrapping channels for one configuration."""
    vertices = length * length
    union_find = HomologyUnionFind(vertices, (length, length))
    for index, pair in enumerate(pairs):
        if (mask >> index) & 1:
            union_find.add_edge(*pair.primal)
    seen: set[int] = set()
    ranks = []
    directions = [0, 0]
    for vertex in range(vertices):
        root, _, _ = union_find.find(vertex)
        if root in seen:
            continue
        seen.add(root)
        component = union_find.component(vertex)
        ranks.append(component.rank)
        if component.direction_0:
            directions[0] = 1
        if component.direction_1:
            directions[1] = 1
    top = max(ranks, default=0)
    return {"wrap_either": int(top > 0), "wrap_cross": int(top == 2),
            "wrap_direction_0": directions[0]}


def enumerate_torus(length: int) -> dict[str, Any]:
    """Every configuration of an ``L x L`` square-bond torus, with its scores.

    Doubled integers are carried rather than halves, so that ``T``, ``B_even``
    and every covariance stay in exact integer arithmetic until the final
    division.  ``twice_q_score`` is ``2T``; ``twice_b_even`` is ``2 B_even``.
    """
    _require(length >= 2, "L must be at least 2")
    pairs = square_bond_pairs(length)
    bonds = len(pairs)
    _require(bonds == 2 * length * length, "unexpected bond count")
    rows: list[dict[str, int]] = []
    for mask in range(1 << bonds):
        primal = configuration_statistics(length, mask, pairs)
        dual = configuration_statistics(length, geometric_dual_mask(mask, pairs), pairs)
        twice_q = 2 * primal["components"] + primal["open_edges"]
        twice_q_dual = 2 * dual["components"] + dual["open_edges"]
        betti = primal["components"] + primal["cycle_rank"]
        betti_dual = dual["components"] + dual["cycle_rank"]
        rows.append({
            "mask": mask,
            **primal,
            "dual_components": dual["components"],
            "dual_open_edges": dual["open_edges"],
            "twice_q_score": twice_q,
            "twice_q_score_dual": twice_q_dual,
            "betti_sum": betti,
            "betti_sum_dual": betti_dual,
            "twice_b_even": betti + betti_dual,
            "x_source": primal["span_rank"] - 1,
            **wrapping_flags(length, mask, pairs),
        })
    return {"length": length, "bonds": bonds, "vertices": length * length, "rows": rows}


def exact_identity_gate(table: Mapping[str, Any]) -> dict[str, Any]:
    """``T - T* = X`` and its two algebraic consequences, configuration by configuration.

    The wrong object this would catch is a dual transport that is not the
    duality map -- naive bit-complement, for instance, which does not swap the
    primal and dual wrapping channels and would break this identity everywhere.
    """
    vertices = table["vertices"]
    failures = {"t_minus_t_dual_is_x": 0, "t_is_half_v_plus_half_b": 0,
                "b_minus_b_dual_is_two_x": 0, "span_rank_differs_from_max": 0}
    x_values: set[int] = set()
    for row in table["rows"]:
        if row["twice_q_score"] - row["twice_q_score_dual"] != 2 * row["x_source"]:
            failures["t_minus_t_dual_is_x"] += 1
        if row["twice_q_score"] != vertices + row["betti_sum"]:
            failures["t_is_half_v_plus_half_b"] += 1
        if row["betti_sum"] - row["betti_sum_dual"] != 2 * row["x_source"]:
            failures["b_minus_b_dual_is_two_x"] += 1
        if row["span_rank"] != row["max_component_rank"]:
            failures["span_rank_differs_from_max"] += 1
        x_values.add(row["x_source"])
    return {
        "configurations": len(table["rows"]),
        "failures": failures,
        "x_values_realised": sorted(x_values),
        "what_span_rank_agreement_means": (
            "two disjoint non-contractible cycles on a torus must be parallel, so "
            "the image rank equals the largest component's rank. That is an "
            "argument; this counts the configurations where it fails, and the "
            "count is the evidence"
        ),
    }


def _mean(values: Sequence[int]) -> Fraction:
    return Fraction(sum(values), len(values))


def _covariance(left: Sequence[int], right: Sequence[int]) -> Fraction:
    count = len(left)
    return (Fraction(sum(a * b for a, b in zip(left, right)), count)
            - _mean(left) * _mean(right))


def covariance_split(table: Mapping[str, Any]) -> list[dict[str, Any]]:
    """``Cov(O,T) = (1/2)Cov(O,B_even) + (1/2)Cov(O,X)``, exactly, per observable.

    The wrong number this stops us believing is a Betti piece and a topological
    piece that do not add back up to the score they were split from -- which
    would mean the two channels being compared in scale are not a decomposition
    of anything.
    """
    rows = table["rows"]
    twice_q = [row["twice_q_score"] for row in rows]
    twice_b_even = [row["twice_b_even"] for row in rows]
    x_source = [row["x_source"] for row in rows]
    out = []
    for name in OBSERVABLES:
        values = [row[name] for row in rows]
        # factors of two: Cov(O, T) = Cov(O, 2T)/2, and B_even = (2B_even)/2.
        total = _covariance(values, twice_q) / 2
        betti = _covariance(values, twice_b_even) / 4
        topological = _covariance(values, x_source) / 2
        out.append({
            "observable": name,
            "cov_with_q_score": str(total),
            "betti_even_piece": str(betti),
            "ambient_homology_piece": str(topological),
            "split_is_exact": total == betti + topological,
            "topological_fraction": (float(topological / total) if total != 0 else None),
        })
    return out


def _bond_radius_squared(length: int, index: int) -> int:
    """Squared torus distance of a bond's midpoint from the origin, in half-units.

    Bond ``2*(y*L+x)`` is the horizontal bond at ``(x, y)`` and ``2*(y*L+x)+1``
    the vertical one, so midpoints sit at ``(2x+1, 2y)`` and ``(2x, 2y+1)`` on a
    lattice of period ``2L``.
    """
    cell, orientation = divmod(index, 2)
    y, x = divmod(cell, length)
    midpoint = (2 * x + 1, 2 * y) if orientation == 0 else (2 * x, 2 * y + 1)
    period = 2 * length
    total = 0
    for coordinate in midpoint:
        wrapped = coordinate % period
        total += min(wrapped, period - wrapped) ** 2
    return total


def spatial_filtration_levels(length: int) -> list[list[int]]:
    """Nested edge sets, growing outward from the origin by torus radius."""
    bonds = 2 * length * length
    radii = sorted({_bond_radius_squared(length, index) for index in range(bonds)})
    levels = []
    for radius in radii:
        levels.append([index for index in range(bonds)
                       if _bond_radius_squared(length, index) <= radius])
    _require(levels[-1] == list(range(bonds)), "the filtration must end at every bond")
    for earlier, later in zip(levels, levels[1:]):
        _require(set(earlier) < set(later), "the filtration must be strictly nested")
    return levels


def _conditional_means(values: Sequence[Fraction], masks: Sequence[int],
                       revealed: Sequence[int]) -> list[Fraction]:
    """``E[Y | F_j]`` for every configuration, exactly.

    At ``p = 1/2`` every configuration is equiprobable, so the conditional mean
    is the plain average over the configurations sharing the revealed bits.
    """
    selector = 0
    for index in revealed:
        selector |= 1 << index
    totals: dict[int, list[Fraction]] = {}
    for value, mask in zip(values, masks):
        key = mask & selector
        bucket = totals.setdefault(key, [Fraction(0), Fraction(0)])
        bucket[0] += value
        bucket[1] += 1
    return [totals[mask & selector][0] / totals[mask & selector][1] for mask in masks]


def scale_decomposition(table: Mapping[str, Any],
                        coordinates: Sequence[str]) -> dict[str, Any]:
    """``Gamma_j = E[D_j D_j^T]`` over a nested spatial filtration, exactly.

    #256's construction: ``m_j = E[Y|F_j]``, ``D_j = m_j - m_(j-1)``, and
    ``Cov(Y) = sum_j Gamma_j``.  Doing it on the joint vector rather than on one
    coordinate at a time is what gives the cross entries an unambiguous
    scale-resolved meaning, because the whole matrix is PSD at every level.

    The wrong result this stops us believing is a per-scale cross-covariance that
    does not sum back to the total, which would mean the scale profile and the
    total covariance are describing different objects.
    """
    rows = table["rows"]
    masks = [row["mask"] for row in rows]
    vectors = {name: [Fraction(row[name]) for row in rows] for name in coordinates}
    levels = spatial_filtration_levels(table["length"])
    previous = {name: [_mean([row[name] for row in rows])] * len(rows)
                for name in coordinates}
    size = len(coordinates)
    gammas: list[dict[str, Any]] = []
    running = [[Fraction(0)] * size for _ in range(size)]
    for level, revealed in enumerate(levels):
        current = {name: _conditional_means(vectors[name], masks, revealed)
                   for name in coordinates}
        increments = {name: [a - b for a, b in zip(current[name], previous[name])]
                      for name in coordinates}
        gamma = [[_mean_product(increments[coordinates[i]], increments[coordinates[j]])
                  for j in range(size)] for i in range(size)]
        for i in range(size):
            for j in range(size):
                running[i][j] += gamma[i][j]
        gammas.append({
            "level": level,
            "revealed_bonds": len(revealed),
            "gamma": [[str(value) for value in row] for row in gamma],
            "is_positive_semidefinite": _is_psd(gamma),
        })
        previous = current
    total = [[_covariance([row[coordinates[i]] for row in rows],
                          [row[coordinates[j]] for row in rows])
              for j in range(size)] for i in range(size)]
    return {
        "coordinates": list(coordinates),
        "levels": gammas,
        "telescopes_to_the_total_covariance": running == total,
        "total_covariance": [[str(value) for value in row] for row in total],
    }


def _mean_product(left: Sequence[Fraction], right: Sequence[Fraction]) -> Fraction:
    return Fraction(sum(a * b for a, b in zip(left, right))) / len(left)


def _is_psd(matrix: Sequence[Sequence[Fraction]]) -> bool:
    """Exact PSD test by leading principal minors of every principal submatrix."""
    size = len(matrix)
    for subset in range(1, 1 << size):
        indices = [i for i in range(size) if (subset >> i) & 1]
        sub = [[matrix[i][j] for j in indices] for i in indices]
        if _determinant(sub) < 0:
            return False
    return True


def _determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    size = len(matrix)
    work = [list(row) for row in matrix]
    result = Fraction(1)
    for column in range(size):
        pivot = next((r for r in range(column, size) if work[r][column] != 0), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            result = -result
        result *= work[column][column]
        inverse = Fraction(1) / work[column][column]
        for r in range(column + 1, size):
            factor = work[r][column] * inverse
            if factor:
                for c in range(column, size):
                    work[r][c] -= factor * work[column][c]
    return result


def walsh_cross_spectrum(table: Mapping[str, Any], left: str, right: str) -> dict[str, Any]:
    """Cross-spectrum by Boolean degree, which must sum to the covariance.

    This is the ``rho = 1`` end of #227's noise semigroup, where the operator is
    the identity.  The wrong number it stops us believing is a degree profile
    whose total is not the covariance it claims to resolve.
    """
    bonds = table["bonds"]
    count = 1 << bonds
    rows = table["rows"]
    order = {row["mask"]: row for row in rows}
    _require(len(order) == count, "enumeration is not one row per configuration")

    def transform(name: str) -> list[int]:
        """Unnormalised fast Walsh-Hadamard transform, in exact integers.

        Every observable here is integer valued, so the butterfly stays integral
        and the ``1/count`` normalisation is applied once at the end.  Carrying
        ``Fraction`` through 18 passes over 262144 entries is what made this
        unaffordable at L=3.
        """
        values = [order[mask][name] for mask in range(count)]
        step = 1
        while step < count:
            for start in range(0, count, step * 2):
                for offset in range(start, start + step):
                    a, b = values[offset], values[offset + step]
                    values[offset] = a + b
                    values[offset + step] = a - b
            step *= 2
        return values

    left_hat, right_hat = transform(left), transform(right)
    scale = Fraction(1, count * count)
    raw: dict[int, int] = {}
    for subset in range(1, count):
        degree = bin(subset).count("1")
        raw[degree] = raw.get(degree, 0) + left_hat[subset] * right_hat[subset]
    by_degree = {degree: value * scale for degree, value in raw.items()}
    total = sum(by_degree.values(), Fraction(0))
    direct = _covariance([order[mask][left] for mask in range(count)],
                         [order[mask][right] for mask in range(count)])
    return {
        "left": left, "right": right,
        "by_degree": {str(degree): str(value) for degree, value in sorted(by_degree.items())},
        "sum_over_degrees": str(total),
        "direct_covariance": str(direct),
        "agrees_with_the_direct_covariance": total == direct,
    }


def degree_one_parity(table: Mapping[str, Any]) -> dict[str, Any]:
    """Degree-one Walsh weight of the two typed pieces, and why they differ.

    The dual transport sends a single-bond character to minus the character of
    the crossing bond, ``chi_i -> -chi_sigma(i)``.  So for any ``f``, the
    duality-**even** part keeps ``(f^({i}) - f^({sigma(i)}))/2`` at degree one and
    the duality-**odd** part keeps the sum.  On the square torus every bond is
    equivalent to its crossing partner, so ``f^({i}) = f^({sigma(i)})`` and the
    even piece's degree-one weight vanishes identically while the odd piece's
    doubles.

    Both ingredients are measured rather than argued: the per-bond coefficients
    of ``B_even`` and of ``X``, and whether the crossing-bond symmetry holds.
    The wrong number this stops us believing is a degree-one Betti loading, which
    would break the cleanest separation the two pieces have.
    """
    bonds = table["bonds"]
    count = 1 << bonds
    order = {row["mask"]: row for row in table["rows"]}
    pairs = square_bond_pairs(table["length"])
    primal_index = {pair.primal: index for index, pair in enumerate(pairs)}
    sigma = [primal_index[pair.dual] for pair in pairs]

    def coefficients(name: str) -> list[Fraction]:
        values = [order[mask][name] for mask in range(count)]
        step = 1
        while step < count:
            for start in range(0, count, step * 2):
                for offset in range(start, start + step):
                    a, b = values[offset], values[offset + step]
                    values[offset] = a + b
                    values[offset + step] = a - b
            step *= 2
        return [Fraction(values[1 << index], count) for index in range(bonds)]

    out: dict[str, Any] = {
        "crossing_map_is_a_bijection": sorted(sigma) == list(range(bonds)),
        "pieces": {},
    }
    for name in ("twice_b_even", "x_source"):
        weights = coefficients(name)
        out["pieces"][name] = {
            "all_degree_one_coefficients_vanish": all(value == 0 for value in weights),
            "distinct_values": sorted({str(value) for value in weights}),
            "invariant_under_the_crossing_map":
                all(weights[index] == weights[sigma[index]] for index in range(bonds)),
        }
    out["mechanism"] = (
        "chi_i -> -chi_sigma(i) under the dual transport, so the duality-even "
        "piece keeps the antisymmetric degree-one combination and the odd piece "
        "the symmetric one. Verified exhaustively at L=2 and L=3; not proved here "
        "for general L"
    )
    return out


def assemble(lengths: Sequence[int] = (2, 3)) -> dict[str, Any]:
    payload: dict[str, Any] = {"schema": SCHEMA, "issue": ISSUE, "tori": {}}
    for length in lengths:
        table = enumerate_torus(length)
        coordinates = ("wrap_either", "twice_b_even", "x_source")
        entry = {
            "bonds": table["bonds"],
            "vertices": table["vertices"],
            "configurations": len(table["rows"]),
            "identity_gate": exact_identity_gate(table),
            "covariance_split": covariance_split(table),
            "scale_decomposition": scale_decomposition(table, coordinates),
            "degree_one_parity": degree_one_parity(table),
        }
        if length <= 3:
            entry["walsh_cross_spectrum"] = [
                walsh_cross_spectrum(table, "wrap_either", "twice_b_even"),
                walsh_cross_spectrum(table, "wrap_either", "x_source"),
            ]
        payload["tori"][str(length)] = entry
    payload["not_established"] = [
        "that any lattice response loads on one piece rather than the other. "
        "This gate establishes that the two pieces are a decomposition, nothing more",
        "that B_even is a local energy field. It is a duality-even measure tangent",
        "that a covariance with X demonstrates a topological defect theory",
        "anything about square-site matching-odd H4, which has no canonical lift to "
        "a square-bond torus of side 2 or 3 and is deliberately absent here",
    ]
    return payload


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--lengths", type=int, nargs="+", default=[2, 3])
    args = parser.parse_args(argv)
    payload = assemble(args.lengths)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
