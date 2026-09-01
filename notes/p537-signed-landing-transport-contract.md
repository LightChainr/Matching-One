# P537 complete landing transfer and two-defect decision

## Decision correction

The N9 and N16 landing minors are exact **supporting controls**.  They prove
that raw kernel-reconnection and readout-pivotal coordinates need not be rank
one, but they do not by themselves falsify the ordinary-four-arm
`P4`/root-Schur lemma.  Their missing ingredients are the separately
normalized axis and tilted geometries and the full bilinear root/slope term.

The first corrected interface is `branch_only` at `a3bc80c8`.  The newer
site-fibre branch
[`7051ad83`](https://github.com/LightChainr/Matching-One/blob/7051ad83a788becb0478de6e0c235376821c22f4/notes/p537-finite-landing-transfer-definition.md)
now freezes the same complete normalization and adds an exact obstruction to
using global component IDs as the landing state.  The older N9/N16 values
remain pipeline controls; they are not the final source/thermal minor.

## Exact state-space correction: land at a finite collar

Commit
[`139c7e58`](https://github.com/LightChainr/Matching-One/blob/139c7e5850c5be1e8e17c3d58be806080b1e2b73/notes/p537-global-four-arm-emptiness.md)
proves a finite-torus identity that changes the computation target.  At an
alternating site `z`, let `b_z` say that the two occupied neighbours are in
distinct global occupied components and `w_z` say that the two vacant
neighbours are in distinct global matching components after removing `z`.
Then

```text
Delta_z q = 1-b_z-w_z,
b_z+w_z <= 1.
```

The second line follows from homology-rank monotonicity.  Hence the state with
two globally distinct occupied arms and two globally distinct vacant
separators is empty on every finite torus.  A canonical ordinary-four-arm row
must therefore record four labelled arms on a **finite collar or annulus
before their outer reconnection**.  Global component IDs make the intended
row empty; local alternation without collar identities is only a relaxed
`near_block`.  The outer attachment and rank transition must remain a
separate state.  Zero-filling an unmatched geometry profile is forbidden.

## Complete fibre functional

Let `g` denote the axis or tilted geometry, with pool weight `1/2`, and put

```text
c_axis = 1/Delta4,
c_tilted = -1/Delta4,
y_g = 2 c_g E_g.
```

Use logit thermal coordinate `t`.  At the pooled root let

```text
S = K-Np,
B = S^2-Np(1-p),
M_t = mean_g Cov(q,S),
Y_t = mean_g Cov(y_g,S),
R = Y_t/M_t,
H_g = y_g-Rq.
```

Write the canonical source as `a=sum_lambda a^lambda`, including its physical
pair normalization, and define

```text
jM^lambda = mean_g Cov(q,a^lambda),
beta_lambda = jM^lambda/M_t.
```

Fix every site except the Bernoulli thermal site `z`.  For `i=X_z in {0,1}`
let

```text
u_i = i-p,
S_i = K_minus-(N-1)p+u_i,
b_i = u_i S_i-p(1-p),
Htilde_i = H_i-E_g H_g,
A_i^lambda = a_i^lambda-E_g a^lambda.
```

The indivisible signed fibre is

```text
Phi_(g,lambda,z)
  = sum_i p^i (1-p)^(1-i)
      Htilde_i {A_i^lambda u_i-beta_lambda b_i}.
```

Equivalently,

```text
Phi = p(1-p) [
        Htilde_mid D_z a^lambda
      + (a_mid^lambda-E_g a^lambda) D_z H_g
      - beta_lambda (K_minus-(N-1)p+1-2p) D_z H_g
      ].
```

The three displayed terms are kernel reconnection, readout pivotality, and
the root/slope Schur allocation.  They are one functional.  Taking their
absolute values separately or treating them as independent evidence changes
the question.

## Landing matrix and finite stop rule

A physical landing label must retain the source-port partitions before and
after the flip, the thermal four-arm partition on the finite collar, the
outer attachment, the two rank values, and an off-port extra-contact flag.
Anonymous component names are quotiented only after transporting the physical
C4 action; serialization labels are not a substitute for that action.

For the same label set in both geometries:

1. form the complete fibre sums with global `p`, means, `R`, and
   `beta_lambda`;
2. take the simultaneous C4 orbit average;
3. combine the separately normalized axis and tilted matrices with the P4
   signs already carried by `y_g`;
4. only then calculate all value and first-thermal-jet minors.

A nonzero final minor falsifies the finite exact pure-thermal factorization.
All zero value minors at one point are insufficient: their first thermal
derivatives must also vanish, or the cleared polynomial minors must vanish on
the declared thermal neighbourhood.  A finite nonzero minor does not by
itself exclude an asymptotic rank-one law; its rate must still be compared
with `M_t/A_N`.

The provisional N25 six-block clean-two-bridge result at
[`ec3941b0`](https://github.com/LightChainr/Matching-One/blob/ec3941b03b2694e827db1cba34766a82e6146a5a/experiments/p537-landing-matrix-preflight-20260901/REPORT.md)
has nonzero projected minors and mixed response under its explicit port-level
contract.  It kills that contract's rank-one claim.  It does not yet decide a
finite-collar ordinary label, because the saved Bell partition does not record
off-port branching, collar arm identities or a transition-resolved thermal
landing.

The same exact N25 population has a sharper but still provisional result at
[`139c7e58`](https://github.com/LightChainr/Matching-One/blob/139c7e5850c5be1e8e17c3d58be806080b1e2b73/experiments/p537-cyclic-bridge-jet-20260901/REPORT.md):
the `(fixed-M response, d/dM response)` Wronskian of `clean_same` and
`clean_reversed` is `+3.475061476262754e-12` and excludes zero exactly.  Thus
cyclic bridge order is a second thermal coordinate within that declared port
contract.  It is not yet a canonical site-flip result, and the signed
same-minus-reversed contrast is not a physical reflection-parity sector until
a dihedral action is frozen.

## Two-scale decomposition

Use one C4-invariant quotient distance in each geometry,

```text
r = d_g(x,y),
s = d_g(z,{x,y}).
```

With dyadic indicators `2^j<=r<2^(j+1)` and
`2^k<=s<2^(k+1)`, define the ordinary block

```text
F_ord[N,j,k]
  = mean_g sum_(lambda,z) E_(g,-z)[
      I_j(r) I_k(s) Pi_C4(1_ordinary Phi_(g,lambda,z))].
```

Endpoint, ordinary, and extra-contact labels form a frozen disjoint
partition, so that

```text
T_t,N = F_endpoint,N + sum_(j,k) F_ord[N,j,k] + R_extra,N.
```

The three scale regimes are `s << r`, `s comparable to r`, and `s >> r`.
The naive comparable-scale absolute account is
`R^4 pi4(R)^3`; even the triangular diagnostic `alpha4=5/4` makes it grow as
`R^(1/4+o(1))`.  It cannot close original U.

## The next mechanism-changing test

For every conserved outer four-arm landing sector `lambda`, construct the
complete Schur-P4 weight matrix `W_lambda` whose slow indices are endpoint
landing type and pivot/rank transition.  Test the two conditional margins

```text
W_lambda rho_lambda = 0,
pi_lambda^T W_lambda = 0.
```

These are stronger than a zero total sum.  They ask whether the signed weight
is a genuine interaction after separately conditioning on each slow variable.

- If both margins vanish sector by sector, a nonzero remainder must pay two
  landing/coupling defects.  This opens a robust route through
  `R^4 pi4(R) pi5(R)^2` or the corresponding two-six-arm bound, followed by a
  scale-flux/coboundary telescoping argument.
- If only one margin vanishes, only one arm upgrade is available and the
  current qualitative square-lattice exponents do not certify the required
  original-U rate.
- If either margin survives in a noncollision slow sector, freeze that sector
  as the minimal leading signed four-arm carrier.  Stop the automatic
  extra-arm route and turn to a prospective sign/scale prediction for that
  sector.

Thus the immediate P0 is the **complete finite-collar landing matrix plus its
value/thermal jet and two conditional margins**, not another global-component
aggregate, raw distance grid, tiny torus, replay counter, or generic
descriptor.

## Root transport and completion standard

Any exact-`p_c` proof must be transported to the pooled root.  A bounded
near-critical coordinate

```text
|p_N-p_c| N pi4(sqrt(N)) = O(1)
```

controls only comparability of arm probabilities.  It does not transport a
signed cancellation.  Completion therefore requires either uniform landing
and margin estimates throughout the root interval or a derivative bound for
the full functional, including centering, `R(p)`, `beta_lambda(p)`, geometry
weights, and the landing law.  Baseline scales for `M_t`, `R`, `R_t`,
`Y_tt`, and `M_tt` must be stated rather than inferred from raw arms.

No new random stream is required for the finite margin decision.  Existing
archives may provide provisional clean-two-bridge coordinates, but none of
the current CSVs contains all fibre IDs, finite-collar arm labels, outer
attachments and both geometry labels needed to reconstruct the canonical
matrix by joining separate datasets.  The completed global-label aggregate
in the site-fibre worktree has no scored canonical result: its strict global
ordinary row is empty by the theorem above, while its nonempty `near_block`
is only a relaxation.
