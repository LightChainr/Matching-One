# Primitive homology C3 orbit certificate

This closes a bounded exact-algebra slice of parent issue #156.  In the basis
`(1, omega)` with `omega=exp(2*pi*i/3)`, the 60-degree unit acts by

`R = [[1,-1],[1,0]]`.

It has determinant one, preserves `a^2-a*b+b^2`, and satisfies `R^3=-I`.
Consequently it acts with order three on primitive unoriented homology lines.
The certificate enumerates every primitive line under a declared invariant
norm cutoff and verifies that these lines split into disjoint three-cycles.

The three character projectors are implemented through exponents of a formal
cube root of unity, so their Gram matrix is exactly `3I` without floating
complex arithmetic.  Under a 60-degree rotation an even spin has C3 charge
`s/2 mod 3`; therefore H4, H8, and H12 occupy charges 2, 1, and 0.

## Boundary

This is an automorphism and character certificate only.  It supplies no
Pinson/Arguin continuum baseline, production covariance, radial score, or
square-site promotion decision.  Parent issue #156 remains open.
