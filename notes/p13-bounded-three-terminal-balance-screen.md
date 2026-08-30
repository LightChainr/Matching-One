# Bounded three-terminal balance-root screen

This analysis applies the exact terminal-reliability evaluator to the complete
connected three-terminal, one-internal-vertex census frozen on `main` at
`e7d63b4`.  The scoring rule was committed separately before the roots were
computed.

For every candidate, the formal three-terminal self-dual balance equation is

```text
P(all three terminals connected) = P(all three terminals separate).
```

The probabilities are exact integer Bernstein polynomials in one common
independent-bond probability `p`.  Each connected candidate has one crossing
in `(0,1)` because the all-connected event increases from zero to one while
the all-separated event decreases from one to zero.

## Result

All 11 connected canonical orbits were scored.  The primary subset contains
the four candidates whose sole internal vertex has degree at least three.

| primary candidate | edges | primitive balance polynomial, low degree first | root | distance to `0.5927460507896` |
|---|---:|---|---:|---:|
| three-spoke star | 3 | `[-1,0,3,-1]` | `0.6527036446661393` | `0.0599575938765393` |
| star plus one terminal edge | 4 | `[-1,1,3,-2]` | `0.5` | `0.0927460507896` |
| star plus two terminal edges | 5 | `[-1,2,3,-3,-1,1]` | `0.3745431865507393` | `0.2182028642388606` |
| complete graph on the four vertices | 6 | `[-1,3,3,-5,-3,6,-2]` | `0.2928932188134525` | `0.2998528319761475` |

The closest candidate in the broader connected set is the reducible
three-terminal path with internal degree two.  Its primitive equation is

```text
p^2 + p - 1 = 0,
```

so its root is `(sqrt(5)-1)/2 = 0.6180339887498948...`, still
`0.0252879379602948...` from the square-site reference.  It is outside the
primary structural filter and does not supply a new cell mechanism.

At the frozen ranking, the three-spoke star is therefore the primary retained
candidate.  This is the familiar exact three-terminal star mechanism rather
than a new square-site construction.

## Scientific decision

The bounded family is now used, not merely enumerated.  It does not contain a
close or structurally new square-site candidate.  Expanding more scoring or
canonicalization around the same one-internal-vertex homogeneous
independent-bond space has low information value.

A subsequent exact search should change the mechanism: correlated cells,
additional internal vertices with an explicit irreducibility reduction, or an
actual periodically embedded critical-polynomial basis.  Decimal proximity by
itself is not a promotion gate.

## Boundary

These are exact reliability polynomials and exact dyadic root isolations for
the declared finite census.  No candidate has a certified periodic tiling,
dual-cell correspondence, square-site probability-law preservation, critical
manifold, percolation bound, or exact square-site threshold.  Issue #13 remains
open.

