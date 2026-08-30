# P14 W5 joint relative-dual connectivity and periodic gluing

## Result

The missing finite-cell state is now explicit.  The pair of terminal
partitions alone has only 14 realized states, and 11 of them remain ambiguous
for the spherical-dual output.  Replacing that projection by

```text
(closed-rim attachment mask, connectivity partition of F_0,...,F_3)
```

closes the configuration map exactly.  Across all 256 W5 edge configurations
this gives 192 labelled states, or 41 simultaneous `D4` orbits.  Every enriched
state determines:

- the primal four-terminal partition;
- the disk-relative boundary partition;
- the spherical-dual terminal partition obtained after gluing the four outer
  boundary arcs to one outer-face vertex.

The quotient removes only connectivity-invisible redundancy: 176 states have
one edge realization and 16 have five.  The fivefold cells are the four
spanning trees and the full cycle that induce the same connected partition of
the four internal dual-face vertices.

## Explicit periodic cell complex

There is a natural periodic embedding, so the obstruction is no longer merely
local.  On an even checkerboard square torus, put one W5 hub in every black
square, retain every square-grid edge as that cell's rim, and share the corner
terminals between cells.  Each grid edge belongs to exactly one black W5.

The disk-relative duals glue without ambiguity: the four boundary leaves
around each white square are identified with its white-face dual vertex.  At
`L=4` the exact census is

| graph | vertices | edges | degree census |
|---|---:|---:|---|
| checkerboard W5 primal | 24 | 64 | 8 of degree 4, 16 of degree 6 |
| glued disk-relative dual | 40 | 64 | 32 of degree 3, 8 of degree 4 |

The primal-edge/dual-edge correspondence is a 64-edge bijection and
`V-E+F=0` on the torus.  The same formula was checked at `L=6`.  Per black W5
cell the incidence densities are `V=3`, `E=8`, `F=5`.

Thus this periodic construction is a genuine primal/dual pair but not a graph
self-duality: vertex density and degree multisets already disagree.  The local
spherical W5 orbit exchange therefore does not become periodic W5
self-duality after the boundary arcs are glued.

## Scientific consequence

The finite-state bottleneck is solved, while the scalar self-duality route is
closed for the natural checkerboard realization.  The next useful object is a
stochastic comparison or exact local transformation between the two explicit,
nonisomorphic periodic graphs, using the committed edge bijection and the 192
state connectivity quotient.  Another homogeneous balance root cannot answer
that question.

This is an exact finite-cell result plus one explicit periodic embedding.  It
does not prove that every possible W5 embedding fails, and it supplies no
threshold identity, domination theorem, or rigorous bound.

## Reproduction

```bash
python3 scripts/p14_w5_joint_relative_dual.py
uv run --with pytest python -m pytest -q \
  tests/test_p14_w5_joint_relative_dual.py \
  tests/test_p14_w5_terminal_duality.py
```
