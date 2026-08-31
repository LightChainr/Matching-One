# Exact terminal-order disk-embedding certificate

Abstract genus zero is not enough to use a four-terminal graph as a planar
hyperedge generator. This control freezes a cyclic terminal order and asks for
an orientable genus-zero rotation system with all terminals on one face in
that order (up to reversal).

For each connected canonical graph with one internal vertex, the checker
exhausts cyclic neighbor orders. Three terminals have one unoriented cyclic
order class. Four terminals have three:

```text
(0,1,2,3), (0,1,3,2), (0,2,1,3).
```

The certificate reports disk-planarity per graph/order pair rather than
silently treating terminal relabelings as the same boundary geometry. A
four-cycle accepts its boundary order but rejects a crossing order; a
four-spoke tree can realize all three orders. `K4` is a direct control showing
that genus zero need not place all terminals on one face.

For every accepted four-terminal witness, all edge subsets are enumerated.
The crossing connectivity `AC|BD` is absent exactly, as required for a disk
embedding with cyclic boundary order `A,B,C,D`.

## Boundary

This is only a bounded embedding and consistency certificate. It does not
construct the planar dual action, derive self-duality equations or a critical
manifold, prove periodic tiling, rank candidates, estimate thresholds, or
prove percolation bounds. Issue 13 remains open.
