#!/usr/bin/env python3
"""Exact adaptive-support escape from the P250 fixed-toggle no-go.

The morphisms are partial and parameter-free.  D/J transfer a non-anchor site
from one typed essential component into the opposite-colour anchor component,
choosing the unique candidate nearest a landing mark.  Ties are left
undefined, rather than broken by a coordinate convention.  This makes the
rule covariant under square-torus isometries and under colour-complement plus
NN/matching hand exchange.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional

try:
    from integer_period_torus import IntegerHomologyUnionFind, axis_integer_torus
except ModuleNotFoundError:
    from scripts.integer_period_torus import IntegerHomologyUnionFind, axis_integer_torus


DEFAULT_OUTPUT = Path("results/exact-p250-adaptive-support-witness/latest.json")


class AdaptiveSupportOracle:
    """Exact component/support oracle on one axis square torus."""

    def __init__(self, L: int):
        if L < 3:
            raise ValueError("L>=3 is required to avoid collided NN images")
        self.geometry = axis_integer_torus(L)
        self.L = L
        self.full_mask = (1 << self.geometry.n) - 1
        self._component_cache: dict[tuple[int, int, bool, bool], Optional[dict]] = {}
        self._support_cache: dict[
            tuple[int, int, int, int, str, bool], Optional[dict]
        ] = {}
        self._leg_cache: dict[tuple[int, int, bool], int] = {}

    def distance_squared(self, first: int, second: int) -> int:
        x0, y0 = self.geometry.coordinates[first]
        x1, y1 = self.geometry.coordinates[second]
        dx = abs(x0 - x1) % self.L
        dy = abs(y0 - y1) % self.L
        dx = min(dx, self.L - dx)
        dy = min(dy, self.L - dy)
        return dx * dx + dy * dy

    def component(
        self,
        mask: int,
        anchor: int,
        *,
        occupied: bool,
        matching: bool,
    ) -> Optional[dict]:
        key = mask, anchor, occupied, matching
        if key in self._component_cache:
            return self._component_cache[key]
        enabled = [
            bool(mask & (1 << vertex)) == occupied
            for vertex in range(self.geometry.n)
        ]
        if not enabled[anchor]:
            self._component_cache[key] = None
            return None
        edges = self.geometry.matching_edges if matching else self.geometry.primal_edges
        union = IntegerHomologyUnionFind(self.geometry.n, self.geometry.periods)
        adjacency = [[] for _ in range(self.geometry.n)]
        for edge in edges:
            if enabled[edge.i] and enabled[edge.j]:
                union.add_edge(edge.i, edge.j, edge.dx, edge.dy)
                adjacency[edge.i].append(edge.j)
                adjacency[edge.j].append(edge.i)
        vertices = {anchor}
        queue = [anchor]
        for vertex in queue:
            for neighbour in adjacency[vertex]:
                if neighbour not in vertices:
                    vertices.add(neighbour)
                    queue.append(neighbour)
        component = union.component(anchor)
        result = {
            "vertices": tuple(sorted(vertices)),
            "size": len(vertices),
            "rank": component.rank,
            "basis": tuple(tuple(vector) for vector in component.basis),
            "colour": "occupied" if occupied else "vacant",
            "graph": "matching" if matching else "NN",
        }
        self._component_cache[key] = result
        return result

    def support(
        self,
        mask: int,
        anchor_D: int,
        anchor_J: int,
        landing: int,
        operation: str,
        *,
        matching_hand: bool,
    ) -> Optional[dict]:
        key = mask, anchor_D, anchor_J, landing, operation, matching_hand
        if key in self._support_cache:
            return self._support_cache[key]
        if operation not in ("D", "J"):
            raise ValueError("operation must be D or J")
        if not (mask & (1 << anchor_D)) or mask & (1 << anchor_J):
            self._support_cache[key] = None
            return None
        occupied = operation == "D"
        matching = matching_hand if occupied else not matching_hand
        source_anchor = anchor_D if occupied else anchor_J
        target_anchor = anchor_J if occupied else anchor_D
        component = self.component(
            mask, source_anchor, occupied=occupied, matching=matching
        )
        if component is None or component["rank"] < 1:
            self._support_cache[key] = None
            return None
        # A physical detach/join transfers the flipped site into the opposite
        # anchor component.  The marked anchors themselves are protected so
        # both ordered branches keep the same typed state contract.
        transferable = []
        for vertex in component["vertices"]:
            if vertex in (anchor_D, anchor_J):
                continue
            bit = 1 << vertex
            changed = mask & ~bit if operation == "D" else mask | bit
            target = self.component(
                changed,
                target_anchor,
                occupied=not occupied,
                matching=not matching,
            )
            if target is not None and vertex in target["vertices"]:
                transferable.append((vertex, target))
        if not transferable:
            self._support_cache[key] = None
            return None
        distances = {
            vertex: self.distance_squared(vertex, landing)
            for vertex, _ in transferable
        }
        minimum = min(distances.values())
        minimizers = tuple(
            vertex for vertex, _ in transferable if distances[vertex] == minimum
        )
        # Undefined ties are the covariant tie-break.  No coordinate order leaks in.
        if len(minimizers) != 1:
            self._support_cache[key] = None
            return None
        selected = minimizers[0]
        target_component = next(
            target for vertex, target in transferable if vertex == selected
        )
        result = {
            "site": selected,
            "site_coordinate": tuple(self.geometry.coordinates[selected]),
            "distance_squared_to_landing": minimum,
            "nearest_minimizer_count": len(minimizers),
            "source_anchor": source_anchor,
            "target_anchor": target_anchor,
            "landing": landing,
            "operation": operation,
            "source_component": component,
            "target_component_after_flip": target_component,
        }
        self._support_cache[key] = result
        return result

    def apply(
        self,
        mask: int,
        anchor_D: int,
        anchor_J: int,
        landing: int,
        operation: str,
        *,
        matching_hand: bool,
    ) -> Optional[tuple[int, dict]]:
        support = self.support(
            mask,
            anchor_D,
            anchor_J,
            landing,
            operation,
            matching_hand=matching_hand,
        )
        if support is None:
            return None
        bit = 1 << support["site"]
        new_mask = mask & ~bit if operation == "D" else mask | bit
        return new_mask, support

    def leg(self, mask: int, root: int, *, matching_hand: bool) -> int:
        key = mask, root, matching_hand
        if key in self._leg_cache:
            return self._leg_cache[key]
        occupied = bool(mask & (1 << root))
        component = self.component(
            mask,
            root,
            occupied=occupied,
            matching=matching_hand if occupied else not matching_hand,
        )
        if component is None:
            raise AssertionError("the root must belong to its colour component")
        value = int(component["rank"] == 1)
        value = value if occupied else -value
        self._leg_cache[key] = value
        return value

    def rectangle(
        self,
        mask: int,
        anchor_D: int,
        anchor_J: int,
        landing: int,
        *,
        matching_hand: bool,
    ) -> Optional[dict]:
        """Return the common-field four-corner rectangle, or None if partial."""

        if len({anchor_D, anchor_J, landing}) != 3:
            return None
        if not (mask & (1 << anchor_D)) or mask & (1 << anchor_J):
            return None
        first_D = self.apply(
            mask,
            anchor_D,
            anchor_J,
            landing,
            "D",
            matching_hand=matching_hand,
        )
        first_J = self.apply(
            mask,
            anchor_D,
            anchor_J,
            landing,
            "J",
            matching_hand=matching_hand,
        )
        if first_D is None or first_J is None:
            return None
        mask_D, support_D = first_D
        mask_J, support_J = first_J
        second_J = self.apply(
            mask_D,
            anchor_D,
            anchor_J,
            landing,
            "J",
            matching_hand=matching_hand,
        )
        second_D = self.apply(
            mask_J,
            anchor_D,
            anchor_J,
            landing,
            "D",
            matching_hand=matching_hand,
        )
        if second_J is None or second_D is None:
            return None
        mask_DJ, support_J_after_D = second_J
        mask_JD, support_D_after_J = second_D
        responses = {
            "L_D": self.leg(mask_D, landing, matching_hand=matching_hand),
            "L_J": self.leg(mask_J, landing, matching_hand=matching_hand),
            "L_DJ": self.leg(mask_DJ, landing, matching_hand=matching_hand),
            "L_JD": self.leg(mask_JD, landing, matching_hand=matching_hand),
        }
        return {
            "base_mask": mask,
            "anchor_D": anchor_D,
            "anchor_J": anchor_J,
            "landing": landing,
            "matching_hand": matching_hand,
            "mask_D": mask_D,
            "mask_J": mask_J,
            "mask_DJ": mask_DJ,
            "mask_JD": mask_JD,
            "supports": {
                "D0": support_D,
                "J0": support_J,
                "J_after_D": support_J_after_D,
                "D_after_J": support_D_after_J,
            },
            "responses": responses,
            "R_plus": (
                responses["L_D"]
                + responses["L_J"]
                - responses["L_DJ"]
                - responses["L_JD"]
            ),
            "R_minus": responses["L_DJ"] - responses["L_JD"],
        }


def _lean_rectangle(row: dict) -> dict:
    """JSON-ready exact certificate without repeating whole vertex sets."""

    supports = {}
    for name, support in row["supports"].items():
        supports[name] = {
            "site": support["site"],
            "site_coordinate": list(support["site_coordinate"]),
            "distance_squared_to_landing": support["distance_squared_to_landing"],
            "nearest_minimizer_count": support["nearest_minimizer_count"],
            "source_component": {
                "size": support["source_component"]["size"],
                "rank": support["source_component"]["rank"],
                "basis": [
                    list(vector) for vector in support["source_component"]["basis"]
                ],
                "colour": support["source_component"]["colour"],
                "graph": support["source_component"]["graph"],
            },
            "target_component_after_flip": {
                "size": support["target_component_after_flip"]["size"],
                "rank": support["target_component_after_flip"]["rank"],
                "basis": [
                    list(vector)
                    for vector in support["target_component_after_flip"]["basis"]
                ],
                "colour": support["target_component_after_flip"]["colour"],
                "graph": support["target_component_after_flip"]["graph"],
            },
        }
    return {
        key: value
        for key, value in row.items()
        if key not in ("supports",)
    } | {"supports": supports}


def exhaustive_l3() -> dict:
    oracle = AdaptiveSupportOracle(3)
    histogram_plus: dict[int, int] = {}
    histogram_minus: dict[int, int] = {}
    valid = 0
    noncommuting_fields = 0
    typed_pairs = 0
    first_witness = None
    for mask in range(1 << oracle.geometry.n):
        for anchor_D in range(oracle.geometry.n):
            if not (mask & (1 << anchor_D)):
                continue
            for anchor_J in range(oracle.geometry.n):
                if mask & (1 << anchor_J):
                    continue
                for landing in range(oracle.geometry.n):
                    row = oracle.rectangle(
                        mask,
                        anchor_D,
                        anchor_J,
                        landing,
                        matching_hand=False,
                    )
                    if row is None:
                        continue
                    valid += 1
                    histogram_plus[row["R_plus"]] = histogram_plus.get(row["R_plus"], 0) + 1
                    histogram_minus[row["R_minus"]] = histogram_minus.get(row["R_minus"], 0) + 1
                    if row["mask_DJ"] != row["mask_JD"]:
                        noncommuting_fields += 1
                    if row["R_minus"] and first_witness is None:
                        first_witness = _lean_rectangle(row)

                    complement = oracle.full_mask ^ mask
                    dual = oracle.rectangle(
                        complement,
                        anchor_J,
                        anchor_D,
                        landing,
                        matching_hand=True,
                    )
                    if dual is None:
                        raise AssertionError("typed partner of a defined rectangle is undefined")
                    if oracle.full_mask ^ row["mask_DJ"] != dual["mask_JD"]:
                        raise AssertionError("typed involution did not exchange ordered finals")
                    if oracle.full_mask ^ row["mask_JD"] != dual["mask_DJ"]:
                        raise AssertionError("typed involution did not exchange reverse finals")
                    if row["R_minus"] != dual["R_minus"]:
                        raise AssertionError("R_minus lost typed-even covariance")
                    if row["R_plus"] != -dual["R_plus"]:
                        raise AssertionError("R_plus lost typed-odd covariance")
                    typed_pairs += 1
    if first_witness is None:
        raise AssertionError("adaptive support rule did not escape the order no-go")
    return {
        "geometry": "axis L=3",
        "defined_rectangles": valid,
        "noncommuting_final_fields": noncommuting_fields,
        "nonzero_R_minus": valid - histogram_minus.get(0, 0),
        "R_plus_histogram": {str(key): value for key, value in sorted(histogram_plus.items())},
        "R_minus_histogram": {str(key): value for key, value in sorted(histogram_minus.items())},
        "typed_partner_rectangles_checked": typed_pairs,
        "first_exact_witness": first_witness,
    }


def support_involution_gate_l3() -> dict:
    oracle = AdaptiveSupportOracle(3)
    checks = 0
    for mask in range(1 << oracle.geometry.n):
        complement = oracle.full_mask ^ mask
        for anchor_D in range(oracle.geometry.n):
            for anchor_J in range(oracle.geometry.n):
                for landing in range(oracle.geometry.n):
                    for operation, dual_operation in (("D", "J"), ("J", "D")):
                        support = oracle.support(
                            mask,
                            anchor_D,
                            anchor_J,
                            landing,
                            operation,
                            matching_hand=False,
                        )
                        dual = oracle.support(
                            complement,
                            anchor_J,
                            anchor_D,
                            landing,
                            dual_operation,
                            matching_hand=True,
                        )
                        site = None if support is None else support["site"]
                        dual_site = None if dual is None else dual["site"]
                        if site != dual_site:
                            raise AssertionError("adaptive support lost typed involution")
                        checks += 1
    return {
        "identity": "support_D(omega,NN)=support_J(C omega,matching), and conversely",
        "defined_and_undefined_support_pairs_checked": checks,
        "tie_rule": "undefined unless the nearest component site is unique",
        "passed": True,
    }


DIHEDRAL = (
    lambda x, y: (x, y),
    lambda x, y: (-x, y),
    lambda x, y: (x, -y),
    lambda x, y: (-x, -y),
    lambda x, y: (y, x),
    lambda x, y: (-y, x),
    lambda x, y: (y, -x),
    lambda x, y: (-y, -x),
)


def vertical_witness_orbit(L: int) -> dict:
    oracle = AdaptiveSupportOracle(L)
    geometry = oracle.geometry
    mask = sum(1 << geometry.vertex((0, y)) for y in range(L))
    anchor_D = geometry.vertex((0, 0))
    anchor_J = geometry.vertex((1, 1))
    landing = geometry.vertex((0, 1))
    base = oracle.rectangle(
        mask, anchor_D, anchor_J, landing, matching_hand=False
    )
    if base is None or base["R_minus"] == 0:
        raise AssertionError("vertical adaptive-support witness failed")

    covariance_checks = 0
    base_supports = {
        name: support["site"] for name, support in base["supports"].items()
    }
    for transform in DIHEDRAL:
        for tx in range(L):
            for ty in range(L):
                def moved(vertex: int) -> int:
                    x, y = geometry.coordinates[vertex]
                    u, v = transform(x, y)
                    return geometry.vertex((u + tx, v + ty))

                moved_mask = sum(
                    1 << moved(vertex)
                    for vertex in range(geometry.n)
                    if mask & (1 << vertex)
                )
                row = oracle.rectangle(
                    moved_mask,
                    moved(anchor_D),
                    moved(anchor_J),
                    moved(landing),
                    matching_hand=False,
                )
                if row is None or row["R_minus"] != base["R_minus"]:
                    raise AssertionError("witness did not survive a square-torus isometry")
                for name, site in base_supports.items():
                    if row["supports"][name]["site"] != moved(site):
                        raise AssertionError("canonical support did not transform covariantly")
                covariance_checks += 1
    return {
        "L": L,
        "base_occupied_coordinates": [[0, y] for y in range(L)],
        "anchors": {
            "D": [0, 0],
            "J": [1, 1],
            "landing": [0, 1],
        },
        "rectangle": _lean_rectangle(base),
        "translation_dihedral_covariance_checks": covariance_checks,
    }


def build_result() -> dict:
    return {
        "schema": "matching-one/p250-adaptive-support-witness/v1",
        "status": "exact_parameter_free_physical_escape",
        "rule": {
            "D": (
                "delete the unique nearest-to-landing non-anchor site in the occupied "
                "D-anchor essential component that joins the vacant J-anchor component after deletion"
            ),
            "J": (
                "add the unique nearest-to-landing non-anchor site in the vacant J-anchor "
                "essential component that joins the occupied D-anchor component after addition"
            ),
            "essential": "ambient homology rank at least one",
            "anchor_contract": "D/J anchors stay occupied/vacant and are never selected",
            "distance": "minimum squared Euclidean distance on the axis square torus",
            "tie_break": "undefined on a tie; never use a coordinate lexicographic choice",
            "typed_state": "(field, NN/matching hand, D anchor, J anchor, landing)",
            "typed_involution": "(omega,h,aD,aJ,c)->(C omega,1-h,aJ,aD,c)",
        },
        "support_involution_gate": support_involution_gate_l3(),
        "exhaustive_L3": exhaustive_l3(),
        "isometric_witness_orbits": [
            vertical_witness_orbit(3),
            vertical_witness_orbit(4),
        ],
        "frozen_minimal_pilot_interface": {
            "base_fields": [
                "replica_id",
                "geometry_id",
                "matching_hand",
                "anchor_D",
                "anchor_J",
                "landing",
                "base_field_digest",
            ],
            "branch_specific_supports": [
                "D0_site",
                "J0_site",
                "J_after_D_site",
                "D_after_J_site",
            ],
            "support_source_descriptors": [
                "colour",
                "graph",
                "component_id",
                "component_size",
                "ambient_rank",
                "primitive_basis",
                "distance_squared_to_landing",
                "nearest_minimizer_count",
            ],
            "responses": ["L_D", "L_J", "L_DJ", "L_JD", "R_plus", "R_minus"],
            "paired_record": "same field under complement+hand exchange with anchors D/J swapped",
            "production_status": "interface frozen by exact witness; no stochastic production run",
        },
        "claim_boundary": (
            "This proves that one explicit covariant adaptive-support rule escapes the fixed-site "
            "order no-go on finite tori. It does not establish path memory in an existing archive, "
            "select this rule as physically unique, or imply a continuum/Jordan statement."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--stdout", action="store_true")
    args = parser.parse_args()
    payload = build_result()
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    if args.stdout:
        print(text, end="")


if __name__ == "__main__":
    main()
