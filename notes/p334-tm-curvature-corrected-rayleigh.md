# Aggregate TM as a curvature-corrected Rayleigh inequality

Fix a carrier, projective line `ell`, and lower layer `F_k(ell)`. Put
`m=N-k`. Choose a state uniformly from the layer and then an ordered pair of
distinct absent sites. Its upward Boolean two-face has exactly four motifs:

- `D` (coexit): both single insertions have rank two;
- `M` (mixed): exactly one single insertion has rank two;
- `Y` (synergy): both singles remain on `ell`, but their double insertion has
  rank two;
- `F` (flat): both singles and the double insertion remain on `ell`.

Write the corresponding integer counts as `D,M,Y,F` and
`T=D+M+Y+F=A m(m-1)`.

## Exact four-motif identity

If `p` is the one-mark exit probability, symmetry of the ordered pair gives

`p = P(D) + P(M)/2`.

Therefore

`P(D) <= p^2 + P(Y)`

is equivalent to

`4DF <= M^2 + 4Y(T-D)`.

Direct integer reduction against the marked-path definition proves

`M^2 + 4Y(T-D) - 4DF = 4(m-1) TM_margin`.

Thus aggregate TM is exactly a curvature-corrected Rayleigh inequality. The
ordinary Rayleigh hard product `D*F` may be covered in two ways:

1. by two mixed faces (`M^2`);
2. by a concave synergy face paired with any non-coexit face (`Y(T-D)`).

A prospective injection has the exact type

`coexit x flat -> mixed x mixed OR synergy x non-coexit`.

## Neither cover can be deleted

Ordinary Rayleigh without the synergy correction passes 968 rows but fails
16. The first failure is matching `diag(2,5)`, `ell=(1,0)`, `k=4`.

Synergy-only coverage passes 916 rows but fails 68. The first failure is a
synergy-free matching row on `diag(3,3)`, `ell=(0,1)`, `k=4`.

The bounded integer cover cone has nine Pareto-minimal ratio points and four
exact lower-convex-hull points. Its two endpoints include a mixed-deficient,
synergy-rich regime and a synergy-free regime. Hence a one-mechanism proof
cannot span the observed topology.

## Stronger local and spectral statements are false

The corrected inequality is not true displacement by displacement. It fails
3,900 of 51,912 individual site-pair tables. The first counterexample is
matching `diag(2,3)`, `ell=(1,0)`, `k=2`, at the opposite order-two pair:

`(T,D,M,Y,F)=(4,2,0,0,2)`, polynomial `-16`.

Grouping every pair only by quotient order still fails 220 of 3,330 tables;
the first order-two aggregate has `(12,6,0,0,6)` and polynomial `-144`.
Thus inversion, torsion order and other delta-local pairings are too strong;
different displacement classes must exchange mass.

A Fourier sum-of-squares is also impossible. With loops doubled in the signed
supply-minus-demand kernel, the exact contrast `e_0-e_1` has quadratic form
`-48` already on primal `diag(2,3)`, `ell=(1,0)`, `k=2`. A simple negative
two-site contrast occurs in 802 rows.

Alexander complement reverses exit faces into birth faces rather than the
reflected exit table. Only 18 of 492 reflected primal/matching margins are
equal; the first pair is `144` versus `288`. Averaging the carriers therefore
does not prove either individual TM inequality.

## The single remaining topology inequality

All 984 aggregate tables satisfy the corrected inequality, but this bounded
fact is not the general proof. The exact remaining statement is:

`4 D F <= M^2 + 4 Y (T-D)`

after summing every relative displacement for one fixed line and layer.
Proving this four-face injection is aggregate TM. Translation regularity then
turns it into the explicit one-mark Hall injection and every proper Hall cut,
with no further hypotheses.
