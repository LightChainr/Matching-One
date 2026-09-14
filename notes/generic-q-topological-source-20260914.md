# Generic-Q topological chemical potential: exact odd shift plus even homology-shape tangent

Date: 2026-09-14. Continuum homology addendum to draft #773 and the graded-Q program. Uses Arguin's critical FK torus homology partition functions. This is not a square-site generic-Q lattice construction.

## 1. Critical FK homology fixes the odd coordinate exactly

For critical FK Potts clusters on a torus, Arguin's homology relation is

    Z_rank2 = Q Z_rank0.

Therefore `P2=Q P0`. In the canonical rank-simplex coordinate

    b = 1/2 log(P0/P2),

the critical FK family lies exactly on

    b_c(Q) = -1/2 log Q                                   (1)

for every torus modulus. A bounded topological source `sX`, `X=rank-1`, translates `b -> b-s`, hence

    s_Q = -1/2 log Q                                      (2)

is the unique source which rebalances the trivial and cross sectors. This odd shift is algebraic/topological, not a fitted continuum amplitude.

## 2. The nontrivial Q response is source-invariant and matching-even

The second coordinate

    d = log[P1/sqrt(P0P2)] = log(2c)                       (3)

is invariant under `sX`. Thus for `0<Q<4` the critical Q tangent decomposes into

    db/d log Q = -1/2                                    (4)

plus the nontrivial

    dd/d log Q.                                           (5)

After applying `s_Q`, the rank law is `(a,r,a)` and the Fisher metric in `(b,d)` is diagonal. Hence (4) and (5) are exactly orthogonal odd/even tangent components of the **aggregate homology law**.

At Q=1, tau=i, the continuum sums give

    dd/d log Q = -0.209344181878...,

with Fisher squared norms approximately

    odd  = 0.15476314,
    even = 0.01033509,

so the even shape piece is about 6.26% of the three-state continuum tangent norm at this modulus. This percentage is not a microscopic-score decomposition until an explicit source map is supplied.

## 3. Interface to the graded Q->1 lattice program

The exact lattice work in #746 separates an ambient-homology/topological score from a duality-even score. The continuum rank law now supplies a necessary response pattern with the same grading:

    forced odd rank response:   db/d log Q = -1/2,
    nontrivial even response:   dd/d log Q.                (6)

This does not identify the microscopic `X` and `B_even` individually with continuum fields. A clean validation pipeline is instead:

1. push the declared finite Q tangent to `(P0,P1,P2)`;
2. transform to `(b,d)`;
3. verify the odd component tends to the exact `-1/2` after the declared Q normalization;
4. compare the residual even component with (5) at the same torus modulus.

A failure in step 3 is a normalization/homology-dictionary failure before any LCFT interpretation.

## 4. Q-dependent critical source-zero collision aspects

The source-invariant

    c=P1/(2 sqrt(P0P2))

controls the two bounded-rank source zeros. On rectangular tori, `c=1` gives a universal static double zero after centering by `s_Q`.

The larger-than-one collision aspects from the Gaussian homology sums are:

| Q | r_*(Q) |
|---:|---:|
|1|1.78782935267966...|
|2|1.86285712512264...|
|3|1.96472673614292...|
|4|2.27122564054142...|

The reciprocal values are equivalent by modular S. These collisions are closure-polynomial exceptional points, not transfer-operator Jordan transitions.

## 5. Square-modulus homology laws

At `tau=i`:

| Q | P0 | P1 | P2 | c |
|---:|---:|---:|---:|---:|
|1|0.3095263|0.3809474|0.3095263|0.6153717|
|2|0.2264243|0.3207272|0.4528486|0.5008040|
|3|0.1863335|0.2546660|0.5590005|0.3945386|
|4|0.1700687|0.1496566|0.6802747|0.2199944|

Applying the exact source (2) rebalances the extreme sectors without changing c.

## 6. Q=4 is an endpoint, not an ordinary tangent point

The real Coulomb-gas parameterisation uses

    sqrt(Q)=2 cos(pi e0/2),

so Q=4 corresponds to `e0=0`, the boundary of the real family. A symmetric numerical derivative in Q would probe Q>4 and becomes complex; it is therefore invalid for the real critical FK family. The committed control records Q=4 finite values and the collision aspect, but marks the Q derivative/Fisher even fraction as `endpoint_not_differentiated`.

This is also the point where logarithmic/marginal finite-size effects are known to require special care. No smooth Q=4 tangent is assumed here.

## 7. New source-centered generic-Q equation of state

For a future generic-Q near-critical homology law, use

    b_Q = b + 1/2 log Q,                                  (7)

so criticality is always `b_Q=0`. Then Q dependence separates into the exact chemical-potential shift already removed by (7) and the genuine deformation of the even shape curve

    d = D_Q(b_Q;tau).                                     (8)

This is cleaner for Q-tangent/logarithmic analysis than differentiating uncentered probabilities.

## 8. Boundaries

- These are continuum critical FK homology weights, not a microscopic generic-Q square-site model.
- Q=4 is not differentiated through the real-family endpoint.
- A homology Q tangent is not automatically an LCFT logarithmic partner.
- No new threshold, Jordan, or original-U claim is made.
