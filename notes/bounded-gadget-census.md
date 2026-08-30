# Bounded one-internal-vertex gadget census

This certificate exhausts finite undirected simple graphs with either three or four terminals and exactly one internal vertex. It quotients the labeled graph space by the full terminal symmetric group using the canonical encoding already on `main`; the internal vertex remains internal.

| terminals | labeled graphs | canonical orbits | connected carriers | connected with internal degree at least 3 |
|---:|---:|---:|---:|---:|
| 3 | 64 | 20 | 11 | 4 |
| 4 | 1,024 | 90 | 58 | 27 |

Every orbit retains its exact labeled multiplicity. Those multiplicities sum back to 64 and 1,024 respectively, so the quotient census is complete within the declared bound. Edge-count histograms are also frozen for the full, connected, and degree-filtered spaces.

The stricter filter requires a connected carrier and degree at least three at the sole nonterminal vertex. This removes isolated, leaf, and degree-two cases. It is deliberately not called a general series-parallel irreducibility test.

## Boundary

No probability polynomial, planarity or periodic-tiling certificate, general series-parallel reduction, self-duality equation, critical manifold, ranking, optimization, or percolation bound is supplied. Issue #13 remains open.
