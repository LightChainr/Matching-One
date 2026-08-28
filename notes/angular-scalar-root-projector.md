# Same-N angular scalar root projector and an L^-7 route

## Status

This is an analytic consequence of the two-sector expansion plus a conditional interpretation of the next matching-odd scalar. It is not a Monte Carlo result and it is not a new threshold estimate.

## 1. Root-level angular projection

For two primitive Gaussian orientations at the same site count `N`, let

\[
c_i=\cos(4\theta_i),\qquad i=1,2,
\]

and let `p_i^*` be the individual finite matching roots. Suppose the central matching function contains

\[
M_i(p_c)=A_0N^{-25/8}+c_iA_4N^{-13/8}+\cdots,
\]

with a common leading thermal slope

\[
M'_i(p_c)=B N^{3/8}(1+o(1)).
\]

The individual roots satisfy, to first order,

\[
p_i^*-p_c=-\frac{A_0}{B}N^{-7/2}-c_i\frac{A_4}{B}N^{-2}+\cdots.
\]

Define the angular scalar and H4 root projectors

\[
P_0[p^*]
=\frac{c_1p_2^*-c_2p_1^*}{c_1-c_2},
\]

\[
P_4[p^*]
=\frac{p_1^*-p_2^*}{c_1-c_2}.
\]

Then the leading H4 root shift cancels from `P0`, while `P4` isolates it:

\[
P_0[p^*]-p_c=-\frac{A_0}{B}N^{-7/2}+\cdots,
\]

\[
P_4[p^*]=-\frac{A_4}{B}N^{-2}+\cdots.
\]

Thus the same-N angular projection converts the ordinary `L^-4` orientation root bias into a potential

\[
\boxed{P_0[p^*]-p_c\sim N^{-7/2}=L^{-7}}
\]

single-size estimator if the next surviving scalar is the proposed `V_<1,4>` contribution.

## 2. Relation to the historical two-size annihilator

The Mertens--Ziff two-size construction removes the leading `L^-13/4` central matching amplitude by combining sizes. The angular scalar projector removes the same H4 sector at fixed `N` by combining orientations.

They are therefore complementary filters:

- size annihilation projects against the leading radial amplitude;
- same-N angular projection projects against the leading spin-4 harmonic.

If both leave an `L^-25/4` scalar residual, both can generate an `L^-7` root correction after division by `M'~L^(3/4)`.

Agreement between the two routes would be much more informative than fitting the same small-L exponent twice.

## 3. Why this is not an immediate Monte Carlo request

The scalar root shift is extremely small. With unit amplitude,

```text
N=65:  N^(-7/2) ~ 4.5e-7
N=85:  N^(-7/2) ~ 1.8e-7
N=185: N^(-7/2) ~ 1.2e-8
```

whereas the existing 100M P45 orientation-root-gap errors are about `7e-6`. The scalar projector retains common-mode absolute-root noise because its weights sum to one. Therefore brute-force Monte Carlo at larger `N` is an inefficient discovery method for this quantity.

Recommended uses are instead:

1. deterministic/exact finite-torus or transfer-matrix root sequences where absolute-root precision is cheap;
2. a future paired estimator only if a pilot demonstrates exceptional common-mode cancellation in `P0[p*]`;
3. a cross-check of an already-detected `L^-7` scalar mechanism, not the primary way to discover it.

Do not allocate a large GPU campaign to this projector without a direct variance pilot.

## 4. Nonlinear and slope caveats

The cancellation above is first-order. Orientation dependence of the thermal slope, H8/H12 contamination, nonlinear root conversion, and logarithmic partners generate subleading leakage. A practical analysis should preserve the full curves and compute `P0[p*]` inside each jackknife replicate rather than projecting already-averaged roots.

The intrinsic pair center `Mbar_N(p0)=0` is not a substitute for this root projector: that center can absorb scalar shifts by construction. The individual roots are required.

## 5. Interpretation boundary

If `P0[p*]` exhibits an `L^-7` law, it supports an angularly scalar correction after the leading H4 sector is removed. Identifying that correction with `V_<1,4>` additionally requires the conditional matching-parity/OPE assumptions in `notes/v14-scalar-post-l7-mechanism.md` and evidence that its lattice coupling is nonzero.
