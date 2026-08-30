# Complement involution on the bounded gadget quotient

Graph complementation commutes with every allowed terminal relabeling. It therefore descends to the canonical terminal-symmetry quotient, is an involution there, and sends a graph with `k` edges to one with `M-k` edges.

| terminals | canonical orbits | self-complementary orbits | two-orbit complement pairs |
|---:|---:|---:|---:|
| 3 | 20 | 0 | 10 |
| 4 | 90 | 2 | 44 |

The orbit accounting is exact: `0+2*10=20` and `2+2*44=90`. The corresponding edge-count histograms are palindromic. For four terminals, the two self-complementary canonical encodings are frozen in the machine artifact.

The test oracle checks commutation for all 64 labeled three-terminal graphs under all 6 terminal permutations and all 1,024 labeled four-terminal graphs under all 24 permutations.

## Boundary

Complementation does not preserve the connected-carrier filter, so no such claim is made. Graph complement is also not being identified with planar or probability duality. No tiling, critical-manifold, ranking, threshold, or bound conclusion is included. Issue #13 remains open.
