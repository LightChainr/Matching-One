# Exact site-threshold maps and the square-site obstruction

Status: bounded literature comparison for Issue #3. This is not a new threshold
calculation or a claim that all possible exact transformations have been excluded.

## Summary

Known exact two-dimensional site thresholds come from preserving a complete
criticality mechanism, not from a decimal coincidence. Three mechanisms are enough
to expose the obstruction for independent site percolation on the square lattice:

1. a lattice is its own matching lattice, so the matching relation fixes
   `p_c = 1/2`;
2. independent site variables are in one-to-one correspondence with independent
   bond variables on another lattice (the line-graph or covering construction);
3. correlations created by a site-to-bond map are confined inside disjoint
   three-terminal cells whose boundary-connectivity law satisfies an exact
   star-triangle/self-duality condition.

Square-site percolation has none of these properties in its standard formulation.
Its matching lattice is the square lattice with both face diagonals, the direct
site-to-bond image has overlapping correlations, and the natural periodic cell is
four-terminal rather than a self-dual three-terminal hyperedge. Generalized critical
polynomials therefore give a convergent finite-basis sequence, not a basis-independent
exact factor.

## Representative exact site models

The table is deliberately short. It lists mechanisms that are structurally relevant
to the square-site question rather than attempting a catalog of every decorated
lattice.

| Site model | Exact condition and threshold | Why the condition applies |
|---|---|---|
| Triangular lattice (and homogeneous self-matching triangulations such as the union-jack lattice) | `p_c = 1/2` | The site matching lattice is the same lattice. The Sykes--Essam threshold relation `p_c(G) + p_c(hat G) = 1` then closes on one unknown. |
| Kagome lattice | `1 - 3p^2 + p^3 = 0`, so `p_c = 1 - 2 sin(pi/18) = 0.652703644666...` | Kagome is the line graph (covering lattice) of the honeycomb lattice. Each honeycomb bond becomes exactly one Kagome site, so independent bond occupation is carried to independent site occupation without changing connectivity. |
| `(3,12^2)` lattice | `p_c^2 = 1 - 2 sin(pi/18)`, so `p_c = 0.807900764120...` | It is the line graph of the doubled-bond honeycomb lattice; the same one-variable-per-edge correspondence applies. |
| Martini lattice | `1 - 3p^3 + p^4 = 0`, so `p_c = 0.764826485929...` | The site-to-bond correlations are encoded inside separated three-terminal cells. Matching their complete boundary-connectivity probabilities through the correlated star-triangle condition gives the exact polynomial. |
| Two martini descendants in Scullard's Figures 8 and 9 | `1-p-p^2=0`, giving `(sqrt(5)-1)/2`, and `1-2p^2=0`, giving `1/sqrt(2)` | Fixing selected site probabilities in the same inhomogeneous martini critical surface produces two other lattices while preserving the solved three-terminal connectivity law. |

The line-graph equivalence is stronger than the informal phrase “turn sites into
bonds.” For a graph `H`, vertices of `L(H)` are edges of `H`, and two vertices of
`L(H)` are adjacent exactly when the corresponding edges of `H` share an endpoint.
Thus an iid Bernoulli variable on every edge of `H` becomes one iid Bernoulli variable
on every vertex of `L(H)`; occupied components correspond exactly.

## Why the same moves do not solve square-site percolation

### The matching identity has two different lattices

For nearest-neighbor square-site percolation, the matching lattice adds both diagonals
of each square face. Consequently,

```text
p_c(square site) + p_c(square matching site) = 1
```

relates two unknown thresholds. It does not reduce to `2 p_c = 1`. This is why the
exact finite matching identities used elsewhere in this repository constrain the
square-site problem without solving its infinite-volume threshold.

### The direct site-to-bond image is correlated

If a bond is declared open precisely when both endpoint sites are occupied, then for
site density `p`

```text
P(edge open) = p^2,
P(two adjacent edges open) = p^3,
```

whereas two independent bonds of marginal `p^2` would have joint probability `p^4`.
Every square-lattice site participates in four edges, so these correlations overlap
from cell to cell. Substituting `p^2` into an independent-bond critical polynomial
therefore changes the probability law and is not an exact site-to-bond reduction.

The martini construction avoids this failure in a specific way: it keeps the required
joint probabilities and arranges the correlated bonds in separated three-terminal
units. Its exactness does not justify discarding correlations on the square lattice.

### The exact three-terminal criterion does not close

For a self-dual 3-uniform hypergraph, criticality can be expressed by equality of the
three boundary vertices being all connected and all disconnected. This is the exact
initial condition reproduced by the critical-polynomial construction on the solvable
class. The square-site model instead has a four-regular, non-self-dual hypergraph
description. Four terminals have additional connectivity partitions; a single
all-connected-versus-none-connected equality does not determine the model.

## What critical-polynomial accuracy does and does not mean

For a finite periodic basis `B`, the generalized critical polynomial compares the
probabilities of two-dimensional and zero-dimensional connectivity. In exactly
solvable cases, the root is independent of basis size (equivalently, a small exact
factor persists). In unsolved cases, the root depends on `B` and is an approximation.

Two published comparisons show both the power and the limitation of the method:

| Problem and construction | Finite calculation | Comparison | Conclusion |
|---|---:|---:|---|
| Kagome **bond** percolation, 36-bond generalized critical polynomial | `0.52440572...` | The contemporary numerical value was `0.52440503(5)`, a difference of about `6.9e-7`. | Very close is not exact; the smallest-basis polynomial was already a known but false conjecture. |
| Square **site** percolation, `6 x 6` rectangular basis | `0.592395070817704` | Modern estimates are near `0.59274605079`, a difference of about `3.51e-4`. | The 2012 finite site polynomials converged much more slowly than comparable bond examples. |

The later semi-infinite-cylinder eigenvalue formulation made the same method family
far more accurate: the 2015 extrapolation reported `0.59274605079210(2)`. The
repository's provenance ledger records a 2024 corrected estimate
`0.5927460507896(1)` from a longer sequence, differing in the twelfth decimal place.
That history is evidence for a powerful estimator and against treating a long decimal
or finite-size eigenvalue equality as an exact critical manifold.

The operational diagnostic is therefore structural:

```text
exact candidate:      the same root/factor survives every compatible basis size;
square-site sequence: the finite roots drift and require extrapolation.
```

## Conclusion and boundary

The square-site obstruction is not merely that its threshold lacks a known closed
form. The known exact routes fail at identifiable interfaces:

- matching sends the model to a different lattice;
- the naive site-to-bond transformation does not preserve independence;
- its correlations are not confined to solved three-terminal cells;
- the four-terminal periodic formulation is not self-dual;
- finite critical-polynomial roots are basis-dependent.

This comparison does **not** rule out a new decorated-cell embedding, correlated
critical manifold, or other exact representation. Such a proposal must preserve the
full probability law and predict a basis-independent identity or another structural
consequence; closeness to `0.592746...` alone is insufficient.

## Primary references

- M. F. Sykes and J. W. Essam, “Exact Critical Percolation Probabilities
  for Site and Bond Problems in Two Dimensions,” *J. Math. Phys.* 5,
  1117 (1964), [doi:10.1063/1.1704215](https://doi.org/10.1063/1.1704215).
- C. R. Scullard, “Exact Site Percolation Thresholds Using the Site-to-Bond
  and Star-Triangle Transformations,” *Phys. Rev. E* 73, 016107 (2006),
  [doi:10.1103/PhysRevE.73.016107](https://doi.org/10.1103/PhysRevE.73.016107),
  [arXiv:cond-mat/0507392](https://arxiv.org/abs/cond-mat/0507392).
- C. R. Scullard, “Percolation Critical Polynomial as a Graph Invariant,”
  *Phys. Rev. E* 86, 041131 (2012),
  [doi:10.1103/PhysRevE.86.041131](https://doi.org/10.1103/PhysRevE.86.041131).
- C. R. Scullard and J. L. Jacobsen, “Transfer Matrix Computation of
  Generalised Critical Polynomials in Percolation,” *J. Phys. A* 45,
  494004 (2012), [arXiv:1209.1451](https://arxiv.org/abs/1209.1451).
- J. L. Jacobsen, “Critical Points of Potts and O(N) Models from Eigenvalue
  Identities in Periodic Temperley--Lieb Algebras,” *J. Phys. A* 48,
  454003 (2015), [doi:10.1088/1751-8113/48/45/454003](https://doi.org/10.1088/1751-8113/48/45/454003),
  [arXiv:1507.03027](https://arxiv.org/abs/1507.03027).
- S. Mertens and R. M. Ziff, “Percolation in Finite Matching Lattices,”
  *Phys. Rev. E* 94, 062152 (2016),
  [doi:10.1103/PhysRevE.94.062152](https://doi.org/10.1103/PhysRevE.94.062152),
  [arXiv:1603.07289](https://arxiv.org/abs/1603.07289).
- Y. Yang and S. Zhou, “Comment on ‘Critical points of Potts and O(N)
  models from eigenvalue identities in periodic Temperley--Lieb algebras’,”
  *J. Phys. A* 57, 258001 (2024),
  [doi:10.1088/1751-8121/ad4d2c](https://doi.org/10.1088/1751-8121/ad4d2c).
