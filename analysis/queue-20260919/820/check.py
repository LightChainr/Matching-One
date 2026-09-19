"""Issue 820: finite support and local alphabet controls; no torus rank kernel."""
from itertools import combinations
from fractions import Fraction
from pathlib import Path
import json


def components(n, edges):
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b); adj[b].add(a)
    labels = [-1]*n
    for root in range(n):
        if labels[root] >= 0:
            continue
        labels[root] = root; stack = [root]
        while stack:
            a = stack.pop()
            for b in adj[a]:
                if labels[b] < 0:
                    labels[b] = root; stack.append(b)
    return labels


def leaf_witness(n, edges, terminals):
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b); adj[b].add(a)
    root = min(terminals); seen = {root}; stack = [root]; tree = set()
    while stack:
        a = stack.pop()
        for b in sorted(adj[a]):
            if b not in seen:
                seen.add(b); stack.append(b); tree.add(tuple(sorted((a,b))))
    assert terminals <= seen
    while True:
        degrees = [0]*n
        for a,b in tree:
            degrees[a] += 1; degrees[b] += 1
        leaves = {a for a in range(n) if degrees[a] == 1 and a not in terminals}
        if not leaves:
            break
        tree = {e for e in tree if not (set(e) & leaves)}
    degrees = [0]*n
    for a,b in tree:
        degrees[a] += 1; degrees[b] += 1
    leaf = min(a for a in terminals if degrees[a] == 1)
    remaining = {e for e in tree if leaf not in e}
    labels = components(n, remaining)
    assert len({labels[a] for a in terminals - {leaf}}) == 1
    assert labels[leaf] != labels[next(iter(terminals - {leaf}))]
    return leaf


def matchings(items):
    if not items:
        yield []
        return
    a = items[0]
    for i in range(1, len(items)):
        b = items[i]
        for rest in matchings(items[1:i] + items[i+1:]):
            yield [(a,b)] + rest


def main():
    support = []
    for n in (3,4,5):
        universe = list(combinations(range(n),2)); tested = 0
        for mask in range(1 << len(universe)):
            edges = [e for i,e in enumerate(universe) if mask >> i & 1]
            labels = components(n, edges)
            for k in range(3,n+1):
                for ts in combinations(range(n),k):
                    if len({labels[a] for a in ts}) == 1:
                        leaf_witness(n, edges, set(ts)); tested += 1
        support.append({'vertices':n,'graphs':1 << len(universe),
                        'connected_terminal_cases':tested,'violations':0})
    pairings = list(matchings(tuple(range(8))))
    allowed = [p for p in pairings if all((a-b)%8 in (1,7) for a,b in p)]
    assert len(pairings) == 105 and len(allowed) == 2
    projector = [Fraction((k-1)*(k-2),2) for k in range(4)]
    assert projector == [1,0,0,1]
    out = {'issue':820,'R2_recommendation':'downgrade broad state-family novelty; preserve narrower rank-image question',
           'terminal_support_controls':support,'degree8_pairings':len(pairings),
           'nonzero_checkerboard_pairings':allowed,'three_edge_block_projector_weights':[str(x) for x in projector],
           'exact_named_scalar_specialization_proved':False,
           'scope':'Small graph support witnesses and local alphabets only; no site-torus enumeration.'}
    Path(__file__).with_name('result.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(out,ensure_ascii=False,indent=2))

if __name__ == '__main__':
    main()
