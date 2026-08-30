# A finite typed rewrite system for aggregate TM

An ordered two-site square above a fixed-line base has only five oriented rank patterns: `D=1222`, `M_left=1212`, `M_right=1122`, `Y=1112`, and `F=1111`. Alexander complement reverses the square and gives `0001`, `0101`, `0011`, `0111`, and `1111` respectively.

In the curvature-corrected Rayleigh polynomial, the unique negative product type is `D x F`. There are exactly two positive rewrite reservoirs:

1. `R_M`: two mixed squares (`M x M`);
2. `R_Y`: a synergy square paired with any non-coexit square (`4Y x nonD`).

Fixing lexicographic token labels and priority `R_M` then `R_Y` gives a terminating, machine-verifiable aggregate rewrite system.

## Extreme-ray witnesses and rewrite regimes

All nine bounded Pareto rays have explicit minimal quotient witnesses; four lie on the exact lower convex hull. Six rays close by `R_M` alone. Three mixed-deficient rays use `R_M` then `R_Y`. Across the full atlas, 968 rows are mixed-only and 16 require synergy rescue.
The 16 rescue rows collapse to 4 exact `(T,D,M,Y,F)` signatures. Even in the most demanding signature, only `133/2880` of the available synergy pool is used after exhausting mixed tokens.

## The unique unclosed general critical pair

The bounded rewrite closes all 984 rows with zero unmatched hard tokens. It is not yet a general topology injection: global lexicographic labels do not tell us how to transform the underlying configurations.

The sole unresolved critical pair is `D x F`. A general rule must cross-switch the two ordered missing-site pairs, or pass through the Alexander-dual birth square, to create `M x M` or `Y x nonD` without image collisions. The exact residual is

`K=max(0,4DF-M^2)-4Y(T-D)`.

All audited rows have `K<=0`. The known N=6 displacement counterexample explains why the rule cannot be delta-local: its hard pair has no mixed or synergy cover in the same displacement class.
