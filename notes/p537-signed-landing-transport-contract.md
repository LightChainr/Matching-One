# P537 joint-incidence landing transport and sector decision

## Decision correction

The N9 and N16 landing minors are exact **supporting controls**.  They prove
that raw kernel-reconnection and readout-pivotal coordinates need not be rank
one, but they do not by themselves falsify the ordinary-four-arm
`P4`/root-Schur lemma.  Their missing ingredients are the separately
normalized axis and tilted geometries and the full bilinear root/slope term.

The first corrected interface is `branch_only` at `a3bc80c8`.  PR #544 is now
**closed without merge**, with recorded head `e8e9c7cf`; its scientific value
line includes
[`e21f6220`](https://github.com/LightChainr/Matching-One/blob/e21f6220775d51450d73d4258c3edf2dbdf0a910/results/p537-finite-collar/REPORT.md).
The strict joint refinement is a separate `branch_only` result at `c958b30b`.
PR #545 is also closed without merge at `e1f19e4c`.  These lifecycle labels
are provenance, not judgments on the mathematics and not permission locks.

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

## Exact full-site-flip result on the relaxed near block

The branch head
[`ddbbad64`](https://github.com/LightChainr/Matching-One/blob/ddbbad6462abb9fcebdaae1f7162e287ed3bebed/results/p537-siteflip-aggregate/summary.json)
now completes the full N25 site-flip enumeration at the common root.  Each
geometry contains `8,388,608` off-site backgrounds and `192,937,984` pair
fibres; the axis and tilted producers took `21.264 s` and `17.995 s`.  The
scorer retains global means, per-source `beta`, the source-absent Schur term,
the fixed-z C4 orbit multiplier and the axis/tilted P4 combination.

The strict global ordinary cell has exactly zero fibres, as required by the
finite-collar theorem.  The separately declared relaxed block—local
alternation, no source-port extra contact, no occupied degree branch in the
preferred witness, but both global landing IDs merged—forms a `15 x 12`
P4/Schur matrix.  Its frozen witness is

```text
rows    = rank 0->1 and rank 1->2, arm_mask=5,
          degree_branch=0, both global landing IDs merged
columns = axial2 source absent / present
det     = -4.039642418513639e-18
```

The exact interval excludes zero.  The search stopped at the first certified
nonzero after `3,301` of `6,930` possible minors.  This kills exact
pure-thermal rank one for that **branch-free relaxed near block**.  It is not
the empty strict global-four-arm row and not yet the finite-collar ordinary
matrix; no first thermal jet or two-margin decision for the collar state has
been delivered.

## Exact common-sector value result and current gap

Closed-unmerged PR #544's scientific commit
[`e21f6220`](https://github.com/LightChainr/Matching-One/blob/e21f6220775d51450d73d4258c3edf2dbdf0a910/results/p537-finite-collar/REPORT.md)
lifts the frozen rows and columns to `B_inf(z,1)`, fixes outer
`J_B=J_W=1`, aggregates corner words, and applies the full root/Schur plus
axis/tilted P4 projection.  Its sole declared minor is
`+2.6904188461441777e-14` with an exact positive interval.  This rejects rank
one for the declared coarse signed-mass matrix.

The four midpoints are

```text
S = [[+2.7100449834e-7, +1.6455142058e-7],
     [-1.6978154702e-7, -3.8139820372e-9]].
```

Hence `S 1=(+4.3555592e-7,-1.7359553e-7)` and
`1^T S=(+1.0122295e-7,+1.6073744e-7)`: every coarse value margin survives.
All selected-row local raw `a`, `qa` and `Ea` sufficient sums are exactly
zero.  The nonzero matrix therefore comes from global centering and
componentwise `beta` Schur allocation, not an observed local raw-source
transfer.  Endpoint unavailable/available is also not source off/on.

The legacy producer stores the source Bell key and collar labels separately.
Branch-only
[`c958b30b`](https://github.com/LightChainr/Matching-One/blob/c958b30b/results/p537-finite-collar-joint/REPORT.md)
therefore reruns the frozen two-site counterfactual fibre with a simultaneous
12-port identity.  It records `1,390` sectors and finds `118` exact-nonzero
sector determinants.  All `106/106` sectors with complete pooled rectangles
and all `30/30` sectors with complete two-geometry rectangles are exact
nonzero.  The lexicographically first full witness has
`det=-4.005552609094306e-18`; all sector sufficient statistics exactly
reconstruct the coarse CSV, matrix and determinant.

The determinant decomposition is highly mixed:

```text
sum_s det(L_s)                    0.2101274%
cross-sector bilinear terms      99.7898726%
det(sum_s L_s)                   2.6904188461441777e-14
```

Thus the parent amplitude is not attributable to one sector, but finite rank
two survives after one common x+y+z component identity is frozen.  The
common-sector **value gate is complete**.

This value result does not by itself record transition adjacency or prove
direct source transmission.  That missing edge gate was subsequently
completed by branch-only
[`22f01e33`](https://github.com/LightChainr/Matching-One/blob/22f01e33c2c386e522b17fc781b3adf70e93548e/results/p537-one-defect-diagonal-edge/REPORT.md).
The stronger kernel-changing graph contains `6,846` exact row classes and
`740,950` physical fibres, with total signed mass `-4.9488399165e-6` and
decision `TWO_INDEPENDENT_DEFECT_GAIN_REJECTED`.  First birth contributes
`-5.8210906659e-6`, while completion contributes `+8.7225074943e-7`.
Within the fixed-anchor, NN-thermal-`z`, source-present kernel-changing set,
the no-contact mask is structurally empty; every selected radius-one edge
touches one or both local thermal arms.  The finite carrier candidate is therefore a
contact-local signed four-arm/OPE functional, not an automatically separated
two-defect object.

The saved edge classes do not retain a common stable row key for the
`c958b30b` joint-sector table, nor do they cover complete original-`U`.
Their signed mass is the fixed-anchor, NN-thermal-`z`, source-present,
alternating, kernel/`g16`-changing logit-thermal submass.  It is not complete
`T_t`, `T_N` or `J_N`.

The later exact N25 stage×contact matrix at `df4a64f6` is
`[-2.88380e-6,-2.93729e-6; -5.32257e-6,+6.19482e-6]`, with determinant
`-3.3498535471e-11`.  Commit
[`f46c38c3`](https://github.com/LightChainr/Matching-One/blob/f46c38c3088c1a9f4df8ab0f256b88639f0b34a3/results/p537-contact-stage-n65/REPORT.md)
then passes the prospective held-out gate on one frozen 20M N65 block:
all four signs are `-- -+`, and `Delta_cs=-8.6882161e-14` has 95% CI
`[-1.1459612e-13,-5.9168198e-14]`.  The decision is
`CONTACT_FUSION_COMPLETION_TRANSMITS`.  The cells are pooled-root allocation
terms; full `J_N` additionally applies the common `A_N/M_t` normalization.
Branch-only `95a695c7` adds full covariance and a post-hoc shape diagnostic
on the byte-identical shards, so it is one dependent reanalysis rather than
another evidence block.

Branch-only `bab37f21` gives a post-reveal noncommuting-completion synthesis.
Current branch head
[`f9ba1ff6`](https://github.com/LightChainr/Matching-One/blob/f9ba1ff6/results/p537-full-t-transport/REPORT.md)
uses the same N65 sufficient statistics to reconstruct the complete response:
`J65=-.00162250989+/-.00018553008`, with carrier share `2.551%` of full
`T65` versus `5.892%` at N25.  Pooled-root displacement changes only `.0155%`
of `|J65|`.  The carrier therefore transmits without saturating the original
response.  All these N65 views remain one dependency group.
The later issue-only gauge audit
([comment](https://github.com/LightChainr/Matching-One/issues/537#issuecomment-5490511026))
shows on exact N9 controls that full `T_t` is invariant under a common thermal
coordinate shift while a filtered contact-stage allocation need not be.  Its
N65 sidecar is not Git-integrated.  Preserve the frozen canonical allocation,
but do not promote it to a coordinate-free operator before #275 fixes the
observable and normalizer.

## Reduced-G4 branch has reached its stop

PR #544 closed without merge at recorded head `e8e9c7cf`.  Its same-named
branch continued after closure to current branch-only head `f9ba1ff6`; the
reduced-L8 artifact remains pinned at
[`22f01e33`](https://github.com/LightChainr/Matching-One/blob/22f01e33c2c386e522b17fc781b3adf70e93548e/results/p537-aggregate-wedge-l8-mc/REPORT.md).
The frozen L8 production gives

```text
G4(8) = 0.0006802040610058551
SE    = 0.0000126970016013406
decision = UNRESOLVED_MODEL_BOUNDARY
```

The registered three-way test is unresolved.  The `L^-35/8` pattern is a
post-hoc discovery from the opened L8 result, not a prospective exponent
certificate; the L6 point also shows a non-monotone correction.  No further
reduced-`G4` size is next.  The normalization-free
`Xi_L=N^2 G4(L) sqrt(P_L)` has already been calculated from existing root
data: L4/L5 are `8.00214/7.99820`, L6 is `7.72546+/-0.08469`, and L8 is
`8.03977+/-0.15006`.  It is reduced-interface context and does not displace
#275's observable/normalizer P0.

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

#537 is now P1 theory.  The remaining target is a uniform near-critical
transport statement for the **nonlocal remainder**, not another local contact
descriptor.  Before any larger-size production, #275 must freeze the original
observable class and physical normalizer and supply opposing prediction
vectors.  Only then may a future N85/N130 frozen-cell block be used as a
prospective discriminator.

The all-`z` exact N25 ledger remains optional P1 accounting and producer
regression: it can replay the 22f slice, c958 restriction and full `J25/A25`
and close `F=C+R`, but it cannot establish asymptotic dominance.  A future
N130 partition must retain every root/Schur term and aligned covariance, yet
it is not in the active production queue.  Priority allocates attention; it
is not a lock or permission system.

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

`c958b30b` completed the joint-sector value decision, `22f01e33` completed the
edge decision, `f46c38c3` completed held-out N65 transmission, and `f9ba1ff6`
measured the small carrier share in full N65 response.  The remaining #537 P1
object is the nonlocal remainder/near-critical theorem; #275 owns the P0
observable and identifiability decision.

## Independent local-collar regression and the missing incidence key

Closed-unmerged #545 at recorded head `e1f19e4c`; its scientific control
[`f6b4414a`](https://github.com/LightChainr/Matching-One/blob/f6b4414ace201846f7c164aea0738a9941eb8d75/experiments/p337-landing-minor-exact-20260901/REPORT.md)
provides a complementary exact regression.  It uses the full canonical pair
source and Bernoulli midpoint thermal score, but one square geometry and the
local collar character `ell4=1_axis-1_diagonal`.  The exhaustive L4 matrix has

```text
det = -533831111/140737488355328 != 0,
```

and an explicit family keeps the analogous two-fibre minor nonzero for every
declared `L,R` at `p>=1/2`; its Sturm certificate also excludes a determinant
zero throughout the exact L4 matching-root interval and gives midpoint
`-2.501463041122436e-6`.  This rejects a single-geometry local-`ell4`
readout common-thermal-multiple and shows that the rank `0->1`/`1->2` split is
not a radius-one accident.  It does not contain the axis/tilted P4 pair,
`H_mid D_z a`, the complete fibre-dependent `-beta H B`, all transition
classes, outer-sector conditioning, the first thermal jet or either margin.
Use its exact source, rank and landing outputs as regression fixtures only.

The historical
[joint-incidence audit](https://github.com/LightChainr/Matching-One/issues/537#issuecomment-5489268565)
adds a required field rather than a new observable.  Separately canonicalized
source and thermal-site Bell keys can agree while `Delta_z g` differs; `93`
of `171` N16 groups mix terminal-incidence classes under that lossy schema.
Every collar fibre must therefore retain

```text
canon_global(x_ports + y_ports + z_ports),
terminal incidence and z-component -> source role,
private grade, joint D4/C4 orbit,
finite-collar arm IDs, outer attachment and rank transition.
```

An equivalent transported component map is acceptable.  `c958b30b` has
implemented this identity and completed the value gate.  The later edge table
does not share a stable row key with that sector output, so the complete
transmission ledger must emit both views from one producer rather than
post-joining them.
