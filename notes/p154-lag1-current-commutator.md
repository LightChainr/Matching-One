# P154 lag-one transmission is a gradient–current pairing

**New exact constraint:** after the pooled root and slope terms are included,
the original U response pairs the three signed birth currents with a rank
potential. Its `02` edge weight is exactly the sum of its `01` and `12` edge
weights. Consequently a current circulation is invisible to every endpoint
rank readout. In particular, a purely orientation-antisymmetric direct-birth
current produces **zero total U response**, while its entry/completion readout
responses can be nonzero and equal-and-opposite. No rigid-clock assumption is
used.

This is finite probability algebra for the existing source, not a new source,
descriptor, sampling proposal or continuum identification. Inputs are
[`4daae57e:src/norm4_lagged_source_replay.cpp`](https://github.com/LightChainr/Matching-One/blob/4daae57eef5c945aa050a95cd3d5d5d77582161b/src/norm4_lagged_source_replay.cpp)
and the full-root formula in
[`bde1a51c:notes/p154-prospective-birth-clock-transmission-decision.md`](https://github.com/LightChainr/Matching-One/blob/bde1a51ca95c74448265b670ba0d9a0d87915479/notes/p154-prospective-birth-clock-transmission-decision.md).
The other team's prospective contract/readout split was read as an uncommitted
draft, not adopted as a freeze or authorization.

## 1. Centering before insertion produces a conserved signed current

Fix geometry g and final occupation K. Let A be the uniform `(K−1)`-site
configuration, V a uniform vacant site, and A⁺=A∪{V}. Write r=R(A),
r⁺=R(A⁺), with monotone ranks 0,1,2, and

```text
s(A)=CB(A)+CW(A),   z(A)=s(A)−E[s|K−1,r,g].
```

The rank-stratified configuration tilt followed by uniform insertion has
first-order endpoint response `J_h(K)=E[z h(r⁺)]`. Since `E[z|r]=0`,

```text
J_h(K)=E[z {h(r⁺)−h(r)}].
```

Let `t_ab=E[z 1{r=a,r⁺=b}]`. The complete transition-score rows sum to zero:
`t00=−t01−t02`, `t11=−t12`, `t22=0`. Thus the endpoint mass derivative is

```text
(δP0,δP1,δP2) = (−t01−t02, t01−t12, t02+t12).
```

These are signed derivatives, not nonnegative flows. With `q=r−1` and
`E=(r−1)²`, the already-established birth interface is

```text
e_K=δF1=t01+t02,     x_K=δF2=t12+t02,
Jq(K)=e_K+x_K,      JE(K)=x_K−e_K.
```

It measures *predictable hazard alignment*: if `λ_ab(A)` is the fraction of
vacant sites giving a→b, and `π_a=P(r=a)`, then

```text
t_ab = π_a Cov(s,λ_ab | r=a).
```

The residual event noise `1{r⁺=b}−λ_ab(A)` is orthogonal to the past-measurable
score z. Consequently neither the sign of a current nor a positive U gain
follows from rank monotonicity or positive birth probabilities.

For the uniform insertion operator `P_K f(A)=E[f(A⁺)|A]` and multiplication
by the rank function h on each level, this is the exact typed commutator

```text
J_h(K)=〈z, (P_K M_h^(K) − M_h^(K−1) P_K) 1〉.
```

Equivalently, writing Π_j for conditional expectation given rank at level j,
`J_h=〈s,(I−Π_(K−1))P_K h〉`. It is the rank-coarse-graining defect of insertion
read against this source. It is not a commutator of two physical interventions
and does not establish order memory. Constant hazards within each early-rank
class force all currents to vanish; a null for this single source does not
establish that stronger closure.

## 2. Original U supplies the potential, without a fitted transmission gain

Let `b_K(p)=BinomPMF(K;N,p)`. Thermalization uses `Jq=Σ_K b_K Jq(K)` and
`JE=Σ_K b_K JE(K)`; only b_K is differentiated. For two geometries define
`P4 f=(f_first−f_second)/δcos4`, `bar f=(f_first+f_second)/2` and, at the
same unperturbed pooled root,

```text
L=N^(13/8)/2,  D=bar(q'),  B=P4(E'),  H=P4(E''),  T=bar(q''),
R_U=H−B*T/D,   U=L*B/D.
```

Assume D≠0. Eliminating rootdot and slopedot from the existing full derivative
gives

```text
v = L/D * { P4(JE') − (B/D)*bar(Jq') − (R_U/D)*bar(Jq) }.
```

Set `α_first=1/δcos4`, `α_second=−1/δcos4`, and define fixed baseline weights

```text
η_gK = (L/D) α_g b'_K,
ξ_K  = −L/(2D²) * {B b'_K + R_U b_K},
ψ_gK(r)=η_gK E(r)+ξ_K q(r).
```

Substitution proves the single exact current formula

```text
v = Σ_(g,K) E_g[z_(K−1) {ψ_gK(r⁺)−ψ_gK(r)}]
  = Σ_(g,K) [(ξ_K−η_gK)t01 + 2ξ_K t02 + (ξ_K+η_gK)t12].
```

This includes direct angular response, root relocation and denominator
response. The three edge weights obey `w02=w01+w12`. It is the particular
rank potential selected by the original U; pooled entry/exit magnitudes do
not fix their pairing with its angular, sign-changing thermal weights.
Population identities here do not remove uncertainty in estimated baseline
jets: any numerical use retains the same paired source/baseline omissions.

## 3. Two mechanism constraints beyond entry/completion rectangles

**Exact current blind direction.** At every g,K the replacement

```text
(t01,t02,t12) → (t01−χ, t02+χ, t12−χ)
```

preserves e_K, x_K and every endpoint rank response, hence also v. This is the
one-dimensional kernel of the three-state incidence matrix. It compares
signed currents, not an allowed dynamical rank loop; arbitrary χ need not be
realizable by the fixed microscopic source. Endpoint birth curves and U alone
therefore cannot identify direct versus sequential event-current structure.
The existing separate `01/02/12` records already retain the distinguishing
information; no extra mark is required.

**Direct-current angular screening.** The coefficient of t02 is independent
of geometry. Hence its contribution to v depends only on the *pooled* t02
profile through root/slope feedback, even if its angular contrast is large.
If the full current has `t01=t12=0` and
`t02,first(K)=−t02,second(K)` for every K, then

```text
JE=0 in both geometries,  bar(Jq)=bar(Jq')=0,
rootdot=slopedot=v=0,
v_entry=−(L/D)P4[(Σ_K b_K t02)'],
v_completion=−v_entry.
```

This is an exact route to strong cancelling readout channels without net U
transmission. It neither requires small event currents nor a common birth-time
translation. More generally JE's fixed-p direct-jump cancellation does **not**
permit discarding pooled t02, because that current moves the root and slope.

For a proposed selective mechanism, first remove the *same full-source*
feedback from both readout derivatives. With `a=rootdot` and `dD=slopedot`,

```text
v_entry = −(L/D)P4(e') + [−L*a*P4(F1'')/D + L*P4(F1')*dD/D²],
v_completion = (L/D)P4(x') + [L*a*P4(F2'')/D − L*P4(F2')*dD/D²].
```

An **event-entry-only** claim `t02=t12=0` forces the feedback-subtracted
completion term to zero, not the raw v_completion. Event-completion-only
`t01=t02=0` forces the analogous entry term to zero. Direct-only forces the
two feedback-subtracted terms to sum to zero. These are necessary operational
mechanism constraints; their converses and positive numerical amplitudes do
not follow. The draft weak/entry-selective/completion-selective rectangles
remain nonexhaustive numerical response targets unless such event support or
angular-current hypotheses are declared independently.

Finally, each K re-injects a separately centered source one step before its
own readout. The t_ab(K) are not derivatives of one globally tilted trajectory
law. Do not cumulatively sum them as a birth-distribution perturbation, infer
a telescoping-in-K conservation law, or interpret the current-kernel dimension
as a field count. The exact content is source → predictable birth current →
the specified global-U potential pairing. No replay, sampling, numerical
template fit, server action or test campaign was performed for this note.
