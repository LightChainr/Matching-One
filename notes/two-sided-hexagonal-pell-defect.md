# Exact two-sided hexagonal Pell defect

This certificate closes the bounded arithmetic part of parent issue #103.
It treats the two recurrences around the hexagonal aspect ratio on equal terms:

- `p^2-3q^2=+1`, beginning `(2,1)`;
- `p^2-3q^2=-2`, beginning `(1,1)`.

Both are advanced by multiplication by the fundamental unit `2+sqrt(3)`, or
equivalently `(p,q) -> (2p+3q,p+2q)`.  For the normalized modulus displacement

`delta = p/(2q)-sqrt(3)/2`

and site count `N=2pq`, exact arithmetic gives

`N*delta = p^2-pq*sqrt(3) = eta*p/(p+q*sqrt(3))`,

where `eta` is the Pell residual.  Hence the `+1` family lies above the
hexagonal modulus and tends to `+1/2`, while the `-2` family lies below and
tends to `-1`; the limiting amplitude ratio is exactly `-2`.

The approach direction is exact too:

`(N*delta-eta/2)*(p+q*sqrt(3))^2 = eta^2/2 > 0`.

## Boundary

This is a geometry-only quadratic-arithmetic oracle.  It makes no statement
about E4 amplitudes, root estimates, covariance, or production suitability.
Those physical and statistical parts of parent issue #103 remain open.
