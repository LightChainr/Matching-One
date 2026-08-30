# Bounded canonical gadget planarity certificate

The checker enumerates orientable rotation systems: each vertex receives a cyclic order of its incident darts, faces are the cycles of the induced face permutation, and genus is computed exactly from

```text
V - E + F = 2 - 2g
```

for each connected component. A genus-zero witness certifies planarity; a positive minimum requires exhausting every distinct cyclic neighbor order.

| terminals | all canonical orbits | planar | nonplanar |
|---:|---:|---:|---:|
| 3 | 20 | 20 | 0 |
| 4 | 90 | 89 | 1 |

The only nonplanar four-terminal orbit is `4:5:1111111111`, the complete graph K5. All 7,776 of its rotation systems were checked, and its minimum orientable genus is one. The connected subset is 57 planar plus K5; the internal-degree-at-least-three subset is 26 planar plus K5.

Independent controls certify K5-minus-one-edge as planar and K3,3 as genus one after all 64 rotation systems.

## Boundary

This is a bounded abstract-graph planarity result. It does not construct a periodic tiling, compute reliability, define a planar dual gadget, solve a self-duality/critical-manifold equation, rank candidates, or imply a threshold/bound. Issue #13 remains open.
