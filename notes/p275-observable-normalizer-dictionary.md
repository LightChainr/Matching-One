# Issue #275: original-observable and normalizer dictionary

**Status (2026-09-01).** This note freezes names and maps for the next
existing-data forward-prediction score.  It does not select a continuum field,
authorize new Monte Carlo, or promote issue-only calculations to Git facts.
The priority is attention rather than permission.

The central rule is that the following are different mathematical objects:

```text
local insertion / propagating sector
        -> marked trace or raw numerator
        -> normalized expectation
        -> geometry projection and pooled moving root
        -> original U or its source response.
```

An unnormalized numerator zero, an invariant-boundary sandwich zero, a
one-point charged zero, a trace zero, a normalized-expectation zero and an
original-`U` zero are six distinct claims.  Existing exact controls already
show that they cannot be substituted for one another.

## 1. Original finite-topology state and exact activation map

On an honest periodic square-cell torus let

\[
r=\operatorname{rank}\operatorname{im}
  [H_1(\text{black NN})\to H_1(T^2)]\in\{0,1,2\}.
\]

The original configuration observables are

\[
q(\omega)=r-1\in\{-1,0,1\},\qquad E(\omega)=q(\omega)^2\in\{0,1\}.
\]

Their normalized expectations are

\[
A_{\rm top}=\langle q\rangle=P_2-P_0,
\qquad
E_{\rm top}=\langle E\rangle=P_0+P_2=1-P_1.
\]

`E_top` is an Alexander-even topology coordinate.  It is not
`<q>^2`, a calibrated continuum energy field, or a propagating sector.

For a uniformly ordered site filtration use the exact convention

```text
K1 = K_minus,
K2 = K_plus,
H1 = 1[K >= K1],
H2 = 1[K >= K2],
F1(p) = E[H1],
F2(p) = E[H2].
```

Then, pointwise and in expectation,

\[
q=-1+H_1+H_2,\qquad E=1-H_1+H_2,
\]

\[
P_0=1-F_1,\quad P_1=F_1-F_2,\quad P_2=F_2,
\]

\[
\langle q\rangle=-1+F_1+F_2,\qquad
\langle E\rangle=1-F_1+F_2.
\]

For a same-`N` orientation pair at one common coordinate,

\[
\Delta_4 A_{\rm top}=\Delta_4F_1+\Delta_4F_2,
\qquad
\Delta_4 E_{\rm top}=-\Delta_4F_1+\Delta_4F_2.
\]

All of `q`, `E`, `P_j`, `F_i`, `A_top` and `E_top` are dimensionless.
`K1/K2` have site-count units.  The threshold CDFs are normalized scalar
coordinates, not two transfer-matrix sectors.

Provenance:

- `main_integrated`, `2da58552e261724e3862d214072a19d4a629fab9`,
  `notes/digital-alexander-duality-proof.md`;
- `open_pr` #267, `e5f39c1d14d54fb2097fd047585f532448918f7a`,
  `scripts/analyze_two_activation_h4.py`,
  `results/norm4-two-activation-h4/latest.json` and `latest.md`.

## 2. Six object kinds

### 2.1 `local_charged_insertion`

A finite-support tensor or source with an explicit nontrivial internal index
or projector, before a global trace or neutral contraction.  A concrete
repository example is the C4 average of

\[
\bar K=iP_{[Q-2,2]}i^\dagger
\]

at a vacant four-edge colour vertex.  At `Q=1` its bounded occupation tangent
is written

\[
S_{\rm av}=-N^{-1}\sum_x t_x.
\]

The coupling is site-average normalized: a tensor perturbation `epsilon/N`
has a dimensionless response per unit `epsilon`; an extensive-source
convention multiplies that response by `N`.  The completed finite N25 value is

\[
\partial_\epsilon U=+0.0018155512845251097.
\]

Artifacts:

- `branch_only`, `923f66b979a6b6132875f783106c041ed3c0c1a9`,
  `notes/local-four-port-transmission-result.md`;
- current index in `open_pr` #267:
  `analysis/four_leg_trace_interface.json`;
- unlabelled local controls are `main_integrated`
  `81e27d9ced8b620d2b8f3fa31ee7a40d242f7a06`,
  `notes/c4-local-odd-pivotal-tangent.md`, and
  `89b86f428e67453595e02c99a5456c11f4d8d77b`,
  `results/local-20260829/P100-marked-pivotal-h4/REPORT.md`.

Not equivalent to:

- an unlabelled four-arm or pivotal spatial mark;
- a torus propagating sector;
- a seam trace;
- a normalized one-point expectation;
- a neutral two-insertion composite.

Four arms do not create Potts `[2]` charge without an explicit colour,
projector or defect mark.

### 2.2 `torus_propagating_sector`

A state-space block

\[
\mathcal H_\lambda=P_\lambda\mathcal H,
\qquad T_\lambda=P_\lambda T P_\lambda,
\]

typed by an internal representation and, separately, a homology/deck sector.
It is not yet a scalar observable.  One must additionally specify transfer
direction and length, spatial sector, boundary vectors or invariant trace,
and transfer normalization.

`P_lambda` is dimensionless.  `Tr(P_lambda T^m)` has unnormalized partition
weight and depends on the transfer convention.

The current authoritative correction is issue-only:

- Issue #275 comment
  `https://github.com/LightChainr/Matching-One/issues/275#issuecomment-5490477497`
  gives a Q4 rational control with `Tr(P2 T^2)=72/25` while the
  invariant-boundary sandwich is zero.

The topology-sector foundation is `main_integrated`
`f9bb310a4d4125f8f52ee8fc7f3239ad35fa334b`,
`notes/modular-homology-channel-classification.md`.

Not equivalent to:

- `P_[Q-2,2] 1=0`, which is a boundary-sandwich statement;
- a trace contribution, which also needs spectrum, multiplicity and cycle;
- a local insertion carrying the same representation label;
- a normalized probability.

The historical `ddf41aa` ordinary-boundary selection statement must not be
upgraded to exclusion of the invariant torus trace.

### 2.3 `homology_marked_trace`

Choose a measure lift, homology restriction and seam first.  In the latest
specified site-cluster colour lift,

\[
Z_r(\sigma)=\sum_{A:r(A)=r}
p^{K(A)}(1-p)^{N-K(A)}c_\sigma(A),
\]

where every occupied NN component has one colour and the horizontal seam
applies permutation `sigma`.  At identity,
`c_id(A)=Q^(C_B(A))`.  For this lift only,

\[
Z_0(\sigma)=Z_0(id),\qquad
Z_2(\sigma)=\frac{\operatorname{fix}(\sigma)}{Q}Z_2(id),
\]

\[
Z_q(\sigma)=\left(\frac{Z_2}{Q}-Z_0\right)
+\frac{Z_2}{Q}\chi_{\rm standard}(\sigma),
\]

\[
Z_E(\sigma)=\left(\frac{Z_2}{Q}+Z_0\right)
+\frac{Z_2}{Q}\chi_{\rm standard}(\sigma).
\]

These are unnormalized partition weights.  Character coefficients can be
signed and are not probabilities or population shares.

This site-cluster calculation is `issue_only`: its code, complete histograms
and receipts were supplied as conversation attachment
`matching_one_trace_reset_20260901` but are not Git artifacts.

Related Git artifacts use a different, explicitly stated lift:

- `branch_only`, `977fea9272c780aea19cc47f8d33324c28a1293e`,
  `notes/closed-source-hypergraph-rc-twist-projection.md`;
- `open_pr` #267, `f81f70c8866e9165c2c9a3c0802b270e66133c24`,
  `notes/four-leg-trace-denominator-interface.md`;
- `open_pr` #267, `b5a0ef976e69f9e0ae639af8475928dff9c0149d`,
  `analysis/four_leg_trace_interface.json`;
- adopted numerical result remains `branch_only`,
  `54352b2eefa651ca482ca84837053c792e82c71e`,
  `results/p337-s4-trace-transmission/score/score.json`.

The fixed Q4 trace there is

\[
Z_{22}=(T+3D_2-4C_3)/6.
\]

Not equivalent to:

- the site-cluster lift versus `Kreg`, closed-source or critical FK lifts;
- a local charged insertion;
- the propagating sector before taking a trace;
- a positive rank-one probability;
- a normalized expectation.

### 2.4 `scalar_composite`

An index-contracted neutral product or an algebraic scalar of the original
state.  The original example is `E=q^2`.  A charged positive control is a
neutral contraction `[2] x [2] -> []`.

The fixed local two-insertion contraction gives

\[
G(Q)=\operatorname{Tr}(\bar K^2)
=\frac{Q(Q-3)(3Q^2-9Q+8)}{8(Q-1)(Q-2)}
=\frac1{2(Q-1)}+O(1).
\]

Its units are the product of the constituent insertion units.  A homogeneous
site-average two-insertion coefficient carries `N^-2`; separation,
contraction and connectedness are part of the definition.

Artifact: `open_pr` #267,
`16a6548af2bb090d9cec7b8e9236e31c0199a3f0`,
`notes/local-pair-two-insertion-obstruction.md` and
`results/local-pair-two-insertion/{REPORT.md,latest.json}`.

Not equivalent to:

- the one-insertion response;
- `Cov(t_x,t_y)` of separately closed occupation marks;
- a propagating-sector trace;
- a normalized two-point function before its denominator is stated.

Rescaling each insertion by `sqrt(Q-1)` changes the mechanism: it regularizes
this two-copy limit but kills the previously finite one-insertion response.

### 2.5 `normalized_expectation`

For a completely specified geometry and measure family,

\[
\langle O\rangle_{g,\alpha}
=Z_{O,g,\alpha}/Z_{g,\alpha}.
\]

The numerator and denominator must use the same lift, seam and coupling.
The ratio is dimensionless when both have the same partition units.

The latest issue-only exact counterexample must remain explicit.  At
`L4,Q4,p=1/2`, character-projecting separately normalized twists gives

\[
-\frac{1795621993608}{8788246642854353}
=-0.00020432084653290933\ldots,
\]

whereas projecting the unnormalized `q` numerator first and dividing by the
common identity partition is exactly zero.

The Git-backed Q4 analogue has rank-one `Z22`, zero raw `q/E` insertion
numerators, but a nonzero normalization-mediated original-`U` response

\[
J_{22}=+5.440121494634842\times10^{-6}.
\]

Not equivalent to:

- `project then divide` versus `divide each twist then project`;
- adding unnormalized numerators across geometries;
- a configuration scalar before expectation;
- original `U`, which adds thermal differentiation and root/slope
  normalization.

### 2.6 `derived_root_conditioned_response`

Original `U` and its perturbative responses form a sixth kind.  They are not
ordinary normalized expectations.  Their complete frozen normalization is
given in the next section.

## 3. Four-layer physical normalizer and original U

The overloaded word `normalizer` is split into four fields.

### 3.1 Partition normalizer

Normalize every geometry under its own physical partition function before
pooling or taking a direction contrast.  A twist-dependent denominator is a
different estimand from the physical identity denominator.

### 3.2 Angular normalizer

For stored orientations `f,s`,

\[
\mathcal P_4X=\frac{X_f-X_s}{\Delta_4},\qquad
\Delta_4=\cos4\theta_f-\cos4\theta_s,
\qquad \bar X=(X_f+X_s)/2.
\]

The exact declared Gaussian representatives determine `Delta4`; it is not a
fitted harmonic amplitude.

### 3.3 Pooled root and thermal-slope normalizer

Let

\[
m_g(p)=\langle q\rangle_g,\qquad e_g(p)=\langle E\rangle_g,
\]

\[
M(p)=\bar m(p),\qquad Y(p)=\mathcal P_4e(p).
\]

The pooled root and its slope are

\[
M(p_0)=0,\qquad D=M_p(p_0)\ne0.
\]

The root does not require `m_f=m_s=0`; it permits `m_f=-m_s`.  Recompute it
inside every aligned delete-one transform.

### 3.4 Size prefactor

\[
A_N=N^{13/8}/2.
\]

This is a frozen finite-size convention, not a newly fitted exponent.

### 3.5 Original U

\[
\boxed{
U_N=A_N\frac{Y_p}{D}\bigg|_{p=p_0}
=\frac{N^{13/8}\mathcal P_4[\partial_p\langle E\rangle]}
 {2\,\partial_p\bar{\langle q\rangle}}
 \bigg|_{\bar{\langle q\rangle}=0}.
}
\]

`U` is dimensionless because the two `p` derivatives cancel their coordinate
unit.  It nevertheless carries the prescribed `N^(13/8)=L^(13/4)` scaling
convention and cannot be compared directly with a raw pair/two-point power.

### 3.6 Complete epsilon response

For a perturbation `epsilon`, at fixed `p` define

\[
hM=\partial_\epsilon M,\qquad hY=\partial_\epsilon Y.
\]

Then

\[
\dot p_0=-hM/D
\]

and the complete root/slope response is

\[
\boxed{
\partial_\epsilon U_N=A_N\left[
\frac{hY_p}{D}
-\frac{Y_{pp}hM}{D^2}
-\frac{Y_phM_p}{D^2}
+\frac{Y_pM_{pp}hM}{D^3}
\right].
}
\]

Its unit is inverse coupling.  Every mechanism must first supply `hM`, `hY`
and their `p` jets, then pass through this identical map.  A character
coefficient, local response or raw pair power cannot be compared directly
with `U`.

Provenance:

- `main_integrated`, `2236d36c80c8a466d9317c929bc33e92a7ca9d33`,
  `predictions/norm4_two_generator_transfer_20260829.yaml`;
- `open_pr` #267, `24973db19c1d97575e9768a6dc187aa0cad5a34b`,
  `notes/p40-thermal-clock-source-quotient.md`;
- `open_pr` #267, `3e6157f237242938e1c1b12415bca256b11896b0`,
  `results/norm4-source-thermal/latest.json`.

## 4. Type-rejection rules

A candidate column or comparison fails closed when any of these differ
without an explicit linear map:

1. measure lift or Q path;
2. internal representation or homology/deck label;
3. insertion order or local/bilocal support;
4. unnormalized trace versus normalized expectation;
5. partition denominator;
6. normalization-before-projection order;
7. geometry-specific versus pooled root;
8. `Delta cos(4 theta)` convention;
9. value versus `p`, source or Q derivative;
10. raw threshold count versus probability;
11. size prefactor or source units;
12. dependency group.

The exact same-Q1 path control in
`results/weak-q-path-comparison/REPORT.md` is a warning: the tied and
rank-projected site-RC paths have opposite original-`U` tangents even though
they share the Q1 baseline.  A common endpoint does not identify a derivative.

## 5. Current assets under this dictionary

- `two-activation-h4`: normalized scalar `q/E` activation coordinates, not
  propagating sectors;
- `rho`-child `E_top`: normalized expectation with a complex geometry
  character, not an energy field or colour trace;
- primitive/C3: homology/direction-marked normalized probabilities, not
  automatically internal Potts charge;
- P43/P57 `A/E/C/W`: physical normalized expectations and thermal clocks;
  a shared H4 label is not a shared operator;
- `Z22`: a homology/colour-marked trace entering `U` through normalization,
  not a local insertion;
- local four-port: a local charged source entering `U` through direct
  numerator, root and slope terms, not a seam trace;
- `J_top/J_bulk`: candidate mechanism labels, not raw-observable kinds until
  they have a source, units and an explicit map to `hM/hY`.

## 6. Unique missing object and bounded fallback

The immediate missing object is theoretical rather than computational:

> two live mechanisms must each supply a unit-bearing forward map from their
> correctly typed sector/insertion/trace to the same correlated raw
> coordinates and then through the same physical normalizer.

The current proposed pair is a semisimple two-activation completion plus a
distinct even bulk singlet versus a two-component Jordan/log-pair transport.
Allowed amplitudes remain in both images; score their design rank, nullspace
and covariance-weighted image intersection.

If those images remain rank-degenerate on the existing primitive/C3,
`rho`-child, P43/P57 and thermal/modulus rows, the only permitted missing
physical input is one of:

1. one phase-calibrated second physical rotation; or
2. one explicit modulus relation.

Name the one selected by the rank audit.  Do not add both, another descriptive
observable, or new Monte Carlo before that degeneracy is shown.
