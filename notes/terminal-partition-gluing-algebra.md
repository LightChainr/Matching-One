# Exact terminal-partition gluing algebra

This Phase-1 control turns boundary connectivity into a finite compositional
object without introducing a graph-specific criticality claim.

For four labeled terminals, all 15 restricted-growth-string partitions are
enumerated. Their join is the transitive closure of the two equivalence
relations. The complete `15 x 15` Cayley table is frozen and satisfies, by
exhaustion:

- associativity on all `15^3 = 3,375` triples;
- commutativity and idempotence;
- the discrete partition as identity;
- the all-connected partition as absorber.

A separate two-port tensor glues

```text
left  = (A,B,x,y),
right = (x',y',C,D),
x=x', y=y',
output = (A,B,C,D).
```

Every one of the 225 input partition pairs has exactly one output partition.
An independent tiny-graph realization—one star connector per nontrivial
block, followed by explicit interface edges—reproduces all 225 outputs.

Consequently independent probability signatures compose bilinearly through
the deterministic integer tensor. A rational synthetic control preserves
total mass exactly.

## Boundary

This certificate does not touch the W5 periodic/relative-dual construction in
PR 438. It does not build planar duality or self-duality equations, search
composition words, infer a critical manifold, estimate a threshold, or prove
a percolation bound. Issue 13 remains open.
