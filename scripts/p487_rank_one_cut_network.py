#!/usr/bin/env python3
"""Bounded independent reference for the rank-one cut theorem (Issue #487).

Period columns are encoded by the HNF matrix [[h, shear], [0, height]].
Sites have labels x+h*y. The checker rejects coincident nearest-neighbour
half-edges. It changes no production topology code. A fixed occupied simple
essential cycle is cut using the actual cyclic order E,N,W,S, not a fitted
bipartition of the trigger graph.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import gcd
from typing import Dict, FrozenSet, Iterable, List, Optional, Set, Tuple

Vector = Tuple[int, int]
Edge = Tuple[int, int]
Cycle = Tuple[Tuple[int, ...], Tuple[int, ...], Vector]
DIRECTIONS = ((1, 0), (0, 1), (-1, 0), (0, -1))


class HnfSquareTorus:
    """An embedded square NN reference graph with distinct incident darts."""

    def __init__(self, h: int, shear: int, height: int) -> None:
        if any(type(v) is not int for v in (h, shear, height)):
            raise TypeError("HNF entries must be integers")
        if h <= 0 or height <= 0 or not 0 <= shear < h:
            raise ValueError("require h,height>0 and 0<=shear<h")
        self.h, self.shear, self.height = h, shear, height
        self.n = h * height
        self.coordinates = tuple((x, y) for y in range(height) for x in range(h))
        self.adjacency = tuple(
            tuple(self.vertex(x + dx, y + dy) for dx, dy in DIRECTIONS)
            for x, y in self.coordinates
        )
        if any(u in row or len(set(row)) != 4
               for u, row in enumerate(self.adjacency)):
            raise ValueError("this reference excludes coincident NN darts/loops")

    def vertex(self, x: int, y: int) -> int:
        quotient, ry = divmod(y, self.height)
        return (x - quotient * self.shear) % self.h + self.h * ry

    def winding(self, x: int, y: int) -> Vector:
        if y % self.height:
            raise ValueError("displacement is not a period")
        b = y // self.height
        if (x - self.shear * b) % self.h:
            raise ValueError("displacement is not a period")
        return (x - self.shear * b) // self.h, b

    def validate_sites(self, occupied: Iterable[int]) -> Set[int]:
        sites = set(occupied)
        if any(type(v) is not int or v < 0 or v >= self.n for v in sites):
            raise ValueError("site labels must lie in [0,N)")
        return sites

    def rank_bfs(self, occupied: Iterable[int]) -> int:
        """Ambient image rank by integer lifted spanning-forest potentials."""
        occupied = self.validate_sites(occupied)
        seen: Dict[int, Vector] = {}
        first = None
        for root in sorted(occupied):
            if root in seen:
                continue
            seen[root] = (0, 0)
            queue = [root]
            for u in queue:
                x, y = seen[u]
                for di, v in enumerate(self.adjacency[u]):
                    if v not in occupied:
                        continue
                    dx, dy = DIRECTIONS[di]
                    candidate = x + dx, y + dy
                    if v not in seen:
                        seen[v] = candidate
                        queue.append(v)
                    else:
                        difference = candidate[0] - seen[v][0], candidate[1] - seen[v][1]
                        if difference == (0, 0):
                            continue
                        w = self.winding(*difference)
                        if first is None:
                            first = w
                        elif first[0] * w[1] != first[1] * w[0]:
                            return 2
        return int(first is not None)

    def rank_union_find(self, occupied: Iterable[int]) -> int:
        """Independent east/north edge union-find; no BFS/cut routine is used."""
        occupied = self.validate_sites(occupied)
        parent = list(range(self.n))
        size = [1] * self.n
        potential = [(0, 0)] * self.n

        def find(v: int) -> Tuple[int, Vector]:
            x = y = 0
            while parent[v] != v:
                dx, dy = potential[v]
                x, y = x + dx, y + dy
                v = parent[v]
            return v, (x, y)

        first = None
        for u in sorted(occupied):
            for di in (0, 1):
                v = self.adjacency[u][di]
                if v not in occupied:
                    continue
                ru, pu = find(u)
                rv, pv = find(v)
                dx, dy = DIRECTIONS[di]
                difference = pu[0] + dx - pv[0], pu[1] + dy - pv[1]
                if ru == rv:
                    if difference != (0, 0):
                        if first is None:
                            first = difference
                        elif first[0] * difference[1] != first[1] * difference[0]:
                            return 2
                elif size[ru] >= size[rv]:
                    parent[rv], potential[rv] = ru, difference
                    size[ru] += size[rv]
                else:
                    parent[ru] = rv
                    potential[ru] = -difference[0], -difference[1]
                    size[rv] += size[ru]
        return int(first is not None)

    def essential_cycle(self, occupied: Iterable[int], reverse_search: bool = False) -> Cycle:
        """First essential fundamental cycle; selected from occupied edges only."""
        occupied = self.validate_sites(occupied)
        seen: Dict[int, Vector] = {}
        parent: Dict[int, Optional[int]] = {}
        parent_direction: Dict[int, int] = {}
        for root in sorted(occupied, reverse=reverse_search):
            if root in seen:
                continue
            seen[root], parent[root] = (0, 0), None
            queue = [root]
            for u in queue:
                x, y = seen[u]
                order = (3, 2, 1, 0) if reverse_search else (0, 1, 2, 3)
                for di in order:
                    v = self.adjacency[u][di]
                    if v not in occupied:
                        continue
                    dx, dy = DIRECTIONS[di]
                    candidate = x + dx, y + dy
                    if v not in seen:
                        seen[v], parent[v], parent_direction[v] = candidate, u, di
                        queue.append(v)
                    elif candidate != seen[v]:
                        up, down = [], []
                        vertex = u
                        while vertex is not None:
                            up.append(vertex)
                            vertex = parent[vertex]
                        vertex = v
                        while vertex not in up:
                            down.append(vertex)
                            vertex = parent[vertex]
                        up = up[:up.index(vertex) + 1]
                        vertices = tuple(up + down[::-1])
                        directions = tuple(
                            [(parent_direction[c] + 2) % 4 for c in up[:-1]]
                            + [parent_direction[c] for c in down[::-1]]
                            + [(di + 2) % 4]
                        )
                        wx = sum(DIRECTIONS[d][0] for d in directions)
                        wy = sum(DIRECTIONS[d][1] for d in directions)
                        winding = self.winding(wx, wy)
                        return vertices, directions, winding
        raise ValueError("occupied graph has no essential cycle")


@dataclass
class RootedNetwork:
    """Permanent nodes are negative; original switchable site labels are nonnegative."""

    adjacency: Dict[int, Set[int]]
    vacancies: FrozenSet[int]
    left: int
    right: int

    def connects(self, inserted: Iterable[int] = ()) -> bool:
        inserted = set(inserted)
        if not inserted <= self.vacancies:
            raise ValueError("insertion contains a non-vacant site")
        if self.left == self.right:
            return True
        active = (set(self.adjacency) - self.vacancies) | inserted
        seen, queue = {self.left}, [self.left]
        for u in queue:
            for v in self.adjacency[u]:
                if v not in active or v in seen:
                    continue
                if v == self.right:
                    return True
                seen.add(v)
                queue.append(v)
        return False

    def activate(self, site: int) -> Optional[RootedNetwork]:
        """One uniform-growth update; None is the rank-two cemetery state.

        This keeps the initial cut fixed. It is not a rule that recuts the torus
        at every step. Each remaining random site still has its own label.
        """
        if site not in self.vacancies:
            raise ValueError("cannot activate a non-vacant site")
        permanent = set(self.adjacency) - self.vacancies
        merge = {site} | (self.adjacency[site] & permanent)
        if self.left in merge and self.right in merge:
            return None
        fixed = merge & permanent
        new_root = min(fixed) if fixed else min(permanent) - 1

        def image(v: int) -> int:
            return new_root if v in merge else v

        adjacency: Dict[int, Set[int]] = {}
        for u, neighbors in self.adjacency.items():
            target = image(u)
            adjacency.setdefault(target, set())
            for v in neighbors:
                if image(v) != target:
                    adjacency[target].add(image(v))
        return RootedNetwork(adjacency, self.vacancies - {site},
                             image(self.left), image(self.right))


@dataclass
class CutResult:
    network: RootedNetwork
    cycle: Cycle
    old_components: Dict[int, FrozenSet[int]]
    ports: Dict[int, FrozenSet[int]]
    left_sites: FrozenSet[int]
    right_sites: FrozenSet[int]
    neutral_sites: FrozenSet[int]
    singleton_triggers: FrozenSet[int]

    def contact(self, u: int, v: int) -> bool:
        return (v in self.network.adjacency[u]
                or bool(self.ports[u] & self.ports[v]))

    def minimal_pairs(self) -> Set[Edge]:
        return {tuple(sorted((u, v))) for u in self.left_sites for v in self.right_sites
                if self.contact(u, v)}

    def minimal_triples(self) -> Set[Tuple[int, int, int]]:
        triples = set()
        for v in self.neutral_sites:
            left = [u for u in self.left_sites if self.contact(u, v)]
            right = [w for w in self.right_sites if self.contact(v, w)]
            for u in left:
                for w in right:
                    if not self.contact(u, w):
                        triples.add(tuple(sorted((u, v, w))))
        return triples

    def bicliques(self) -> List[dict]:
        result = []
        for component, old_sites in sorted(self.old_components.items()):
            if component in (self.network.left, self.network.right):
                continue
            left = sorted(v for v in self.left_sites if component in self.ports[v])
            right = sorted(v for v in self.right_sites if component in self.ports[v])
            if left and right:
                result.append(dict(component=component, old_sites=sorted(old_sites),
                                   left=left, right=right, edges=len(left)*len(right)))
        return result


def cut_rank_one(torus: HnfSquareTorus, occupied: Iterable[int],
                 cycle: Optional[Cycle] = None) -> CutResult:
    """Build the fixed-cut component network without evaluating any future pair."""
    occupied = torus.validate_sites(occupied)
    if torus.rank_bfs(occupied) != 1:
        raise ValueError("the cut theorem requires ambient rank one")
    cycle = cycle if cycle is not None else torus.essential_cycle(occupied)
    vertices, directions, winding = cycle
    gamma = set(vertices)
    if len(gamma) != len(vertices) or len(vertices) != len(directions) or len(vertices) < 3:
        raise ValueError("cycle must be simple and have one direction per vertex")
    if not gamma <= occupied or any(type(d) is not int or d not in range(4) for d in directions):
        raise ValueError("cycle is not an occupied embedded edge cycle")
    if any(torus.adjacency[u][d] != vertices[(i+1) % len(vertices)]
           for i, (u, d) in enumerate(zip(vertices, directions))):
        raise ValueError("cycle edge mismatch")
    actual_winding = torus.winding(sum(DIRECTIONS[d][0] for d in directions),
                                  sum(DIRECTIONS[d][1] for d in directions))
    if winding != actual_winding or gcd(abs(winding[0]), abs(winding[1])) != 1:
        raise ValueError("cycle must have the declared primitive nonzero winding")
    position = {u: i for i, u in enumerate(vertices)}
    side = {}
    for i, u in enumerate(vertices):
        outgoing, backward = directions[i], (directions[i-1]+2) % 4
        span = (backward-outgoing) % 4
        for d in range(4):
            if d not in (outgoing, backward):
                side[u, d] = int((d-outgoing) % 4 >= span)

    def endpoint(u: int, d: int) -> int:
        return torus.n + 2*u + side[u, d] if u in gamma else u

    cut_adj: Dict[int, Set[int]] = {}

    def edge(u: int, v: int) -> None:
        cut_adj.setdefault(u, set()).add(v)
        cut_adj.setdefault(v, set()).add(u)

    for u in range(torus.n):
        for copy in ((torus.n+2*u, torus.n+2*u+1) if u in gamma else (u,)):
            cut_adj.setdefault(copy, set())
        for d, v in enumerate(torus.adjacency[u]):
            if u > v:
                continue
            along = u in gamma and d in (
                directions[position[u]], (directions[position[u]-1]+2) % 4)
            if along:
                edge(torus.n+2*u, torus.n+2*v)
                edge(torus.n+2*u+1, torus.n+2*v+1)
            else:
                edge(endpoint(u, d), endpoint(v, (d+2) % 4))

    active = (occupied-gamma) | {torus.n+2*u+i for u in gamma for i in (0, 1)}
    component: Dict[int, int] = {}
    old_components = {}
    for root in sorted(active):
        if root in component:
            continue
        label = -1-len(old_components)
        queue, component[root] = [root], label
        for u in queue:
            for v in sorted(cut_adj[u]):
                if v in active and v not in component:
                    component[v] = label
                    queue.append(v)
        old_components[label] = frozenset(queue)
    left, right = component[torus.n+2*vertices[0]], component[torus.n+2*vertices[0]+1]
    if left == right:
        raise AssertionError("rank-one graph unexpectedly bridges its cut boundaries")
    adjacency: Dict[int, Set[int]] = {}
    for u, neighbors in cut_adj.items():
        a = component.get(u, u)
        adjacency.setdefault(a, set())
        for v in neighbors:
            b = component.get(v, v)
            if a != b:
                adjacency[a].add(b)
    vacancies = frozenset(set(range(torus.n))-occupied)
    ports = {v: frozenset(w for w in adjacency[v] if w < 0) for v in vacancies}
    singletons = frozenset(v for v in vacancies if left in ports[v] and right in ports[v])
    safe = vacancies-singletons
    left_sites = frozenset(v for v in safe if left in ports[v])
    right_sites = frozenset(v for v in safe if right in ports[v])
    return CutResult(RootedNetwork(adjacency, vacancies, left, right), cycle,
                     old_components, ports, left_sites, right_sites,
                     safe-left_sites-right_sites, singletons)


def pair_statistics(pairs: Iterable[Edge]) -> dict:
    """Degree-derived summaries only after geometric edge reconstruction."""
    pairs = set(pairs)
    degrees = Counter(v for pair in pairs for v in pair)
    return dict(edges=len(pairs), nonisolated=len(degrees),
                wedges=sum(d*(d-1)//2 for d in degrees.values()))
