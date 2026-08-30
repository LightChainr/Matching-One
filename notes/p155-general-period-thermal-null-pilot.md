# Issue #155: general-period thermal-null pilot

## Frozen discriminator

This bounded pilot carries the exact microscopic counterterm

```text
O_alpha* = O_local_H4 + (3/64) epsilon_cell
epsilon_cell = [(n_(0,0)-1/2) + (n_(1,0)-1/2)] / 2
```

to the first two parity-preserving Gaussian checkerboard quotients whose
Euclidean `R=8` landing registry is injective:

```text
N260: periods [[16,-2],[2,16]]
N340: periods [[18,-4],[4,18]]
```

The checkerboard graph has every nearest-neighbour edge and the two forward
diagonals from each even `x+y` site.  Both period columns have even coordinate
sum, so parity is well-defined on the quotient.  Quotient labels, winding and
ambient-H1 rank use the already-tested general-period HNF backend.

The exact N10 oracle used the original Chebyshev `R=1` landing and is retained
as a graph/readout regression: the rows `(global, local, epsilon)` respond as
`((15/8,5/4),(-3/64,11/64),(1,0))`.  Production uses Euclidean radii `R=2,4,8`
because Chebyshev `R=8` aliases on N260.  Therefore `3/64` is a frozen UV
counterterm, not a fitted claim that the finite-R thermal response is exactly
zero.

At `p_even=p_odd=1/2`, every configuration is paired with its complement.
The score columns are the exact product-measure derivatives `(S_t,S_lambda)`.
The same seed/counter stream is used at both sizes and all radii.

## Decision rule

Run 20,000 configurations per size, 100 aligned batches, seed
`15583020260830`, counters `[15500000000,15500020000)`.  `R=8` is primary;
`R=2,4` are radial diagnostics.

The R8 response matrix has rows `(global,O_alpha*)` and columns `(t,lambda)`.
A size passes only if all three conditions hold:

```text
abs(O_alpha* lambda response z) >= 3
abs(determinant z) >= 3
condition number <= 50
```

One expansion to 100,000 is allowed only when the 20k wall time is below two
minutes, both R8 central condition numbers are already at most 50, and simple
`sqrt(5)` variance projection takes both z gates above 3 at both sizes.  No
further expansion is permitted.

This experiment can resolve a conditioned second finite-volume response
direction.  It does not identify an RG eigenoperator, exponent or CFT field.

