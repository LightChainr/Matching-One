# Distinct microscopic births can still coalesce: a persistence criterion

Date: 2026-09-29. New finite/probabilistic analysis for the paired-birth route.
This note proves statements about ordered two-birth processes; it does **not**
prove a percolation scaling limit or square-site universality.

## The useful new target

The question is not only whether a single insertion can create both homology
directions. It is whether two **different** insertions can become one birth
after near-critical rescaling. Measure the shrinking neighbourhood of the
diagonal, not just its microscopic zero atom.

For a nonnegative, dimensionless birth gap `G`, define

```text
Z(delta) = E[(1-G/delta)_+],                     delta > 0.       (1)
```

This soft near-diagonal mass has two exact representations: one from the paired
gap histogram, another from one- and two-time rank-one occupation. It provides
a criterion for excluding limit double births without fitting a free exponent.

The relevant prior assets are unmerged #773, commit
`9ef22d1cc531c09dedf9b1195ab6aad8ef79668f`:

- [`notes/rank-birth-flux-conservation-20260914.md`](https://github.com/LightChainr/Matching-One/blob/9ef22d1cc531c09dedf9b1195ab6aad8ef79668f/notes/rank-birth-flux-conservation-20260914.md):
  direct jump flux, static marginals and the missing copula.
- [`notes/canonical-two-birth-scaling-20260914.md`](https://github.com/LightChainr/Matching-One/blob/9ef22d1cc531c09dedf9b1195ab6aad8ef79668f/notes/canonical-two-birth-scaling-20260914.md):
  canonical clock and the proposed insertion scale.
- [#778](https://github.com/LightChainr/Matching-One/issues/778): paired archive
  semantics and exact order-statistic mixture. Its latest September 14
  comments also retain the optional primitive rank-one direction mark.

In particular, `P(D=0)->0` alone does not close the continuity/no-aggregation
step of #773. This note supplies the missing distinction, not a refutation of
the proposed scaling programme.

## 1. A counterexample with identical entire one-time rank laws

For every integer `n>=2`, consider two equally weighted birth pairs in each
of two models:

| Model | First pair `(S,T)` | Second pair `(S,T)` |
|---|---|---|
| A | `(0, 1+1/n)` | `(1, 2)` |
| B | `(0, 2)` | `(1, 1+1/n)` |

Set `R(t)=1{S<=t}+1{T<=t}`. Both are monotone rank processes starting at zero
and ending at two. In **both models, at every finite n**, `S<T`: neither ever
has a direct jump of size two. They have exactly the same birth marginals,
hence the same full one-time law of `R(t)` for every t, and the same mean gap
`E(T-S)=1+1/(2n)`.

But their limit gaps are

```text
A: G_n => 1,
B: G_n => (delta_0 + delta_2)/2.                              (2)
```

Thus even all static rank probabilities **plus an identically zero direct
jump flux** cannot decide whether distinct births aggregate. This is a
counterexample in the class of monotone two-birth processes, not an alternative
percolation model. It can also be put on insertion levels by taking
`J_i=1+n*(S or T)` and `N=2n+1`; all resulting indices are valid integers.
In that discrete version the **complete adjacent-time rank transition tables**
are identical too: the static differences determine entry/exit fluxes once
the common direct-jump flux is zero. Even that one-step kernel is therefore
insufficient for the aggregation question; a Markov reconstruction would
silently discard relevant path dependence.

There is an additional moment trap. A positive gap that equals `1/n` with
probability `1-1/n` and `n` with probability `1/n` converges to zero in
probability while its mean converges to one. Tightness alone does not license
passing occupation means to the limit: uniform integrability is needed.

## 2. The occupation-defect theorem

Let `(S,T)` be any almost surely finite ordered birth pair, `G=T-S>=0`, and
assume `E G<infinity`. Define

```text
q(t)       = P(S<=t<T),
H_delta(t) = P(S<=t, t+delta<T),
A(delta)   = integral_R H_delta(t) dt.
```

Then, exactly,

```text
integral q(t) dt = E G,
A(delta)         = E[(G-delta)_+],
E G - A(delta)   = E[min(G,delta)],
Z(delta)         = 1 - [E G-A(delta)]/delta.                   (3)
```

**Proof.** For each realised interval `[S,T)`, the set of t for which both
t and `t+delta` lie in it has Lebesgue length `(G-delta)_+`. Integrate the
indicators and use Tonelli; the remaining equalities are pointwise identities.
No Markov assumption is used.

The last equality is useful conceptually; numerically evaluate (1) directly
from gaps, rather than subtract two nearly equal large occupation integrals.
Formula (1) is meaningful even when `E G` is infinite, in which case the
subtraction in (3) must not be used.

For `0<c<1`, the following bounds convert the observable into a probability:

```text
P(G=0) <= Z(delta) <= P(G<delta),
P(G<=c*delta) <= Z(delta)/(1-c).                              (4)
```

Hence `Z(delta)` is not the zero atom at a fixed finite resolution. It combines
that atom with short but positive bars. For a fixed distribution,
`lim_{delta down to 0} Z(delta)=P(G=0)` by bounded convergence.

### Noncoalescence criterion for a sequence

Suppose the rescaled ordered pairs `(S_n,T_n)` are tight in `R^2`. The following
are equivalent:

1. Every subsequential weak limit satisfies `P(S=T)=0`.
2. `lim_{epsilon down to 0} limsup_n P(T_n-S_n<=epsilon)=0`.
3. `lim_{delta down to 0} limsup_n Z_n(delta)=0`.

Equivalence of 2 and 3 follows from (4), for example with `c=1/2`.
For 2 implies 1, use the open sets `{t-s<epsilon}` and Portmanteau, then let
epsilon decrease to zero. Conversely, failure of 2 gives a subsequence with
uniform positive mass in shrinking diagonal neighbourhoods; tightness gives
a weakly convergent further subsequence. For every fixed epsilon its mass in
the closed set `{t-s<=epsilon}` is bounded below, so the limit has positive
diagonal mass. This proves 1 implies 2.

Tightness of the gaps alone gives the analogous statement about gap limits;
it does not locate the birth pair in a finite thermal window. Pair tightness
or a separate location control is needed for that stronger statement.

The order of limits is essential. Taking `delta->0` first returns only
`P(G_n=0)` and misses model B in (2). In measurements, take an increasing size
sequence at several **fixed scaled** delta values, then examine smaller delta.
A single finite-size grid cannot prove the double limit.

### Exact discrete version

For insertion births `1<=J1<=J2<=N`, put `D=J2-J1` and
`R_k=1{J1<=k}+1{J2<=k}`. Extend the rank-one indicator by zero outside
`k=0,...,N`. For every nonnegative integer m,

```text
A_m = sum_k P(R_k=1, R_{k+m}=1) = E[(D-m)_+],
A_0 = sum_k P(R_k=1) = E D,
1-(A_0-A_m)/m = E[(1-D/m)_+],                     m>=1.       (5)
```

At microscopic resolution `m=1`, this equals `P(D=0)`. At a mesoscopic
resolution `m approximately delta * insertion_window`, it tests aggregation.
The whole distribution can be recovered by discrete curvature:

```text
P(D=d) = A_{d-1} - 2 A_d + A_{d+1},                d>=1.      (6)
```

Thus the new object is a lagged rank-one overlap, not a larger catalogue of
one-time observables. An existing paired histogram already determines it.

## 3. Exactly what marginal information can establish

Positive mean rank-one occupation is useful, but it does not exclude a partial
diagonal atom; models A and B have the same mean and all the same marginals.
With `m=E G` and `v=E G^2<infinity`, Cauchy--Schwarz does give

```text
P(G>epsilon) >= (m-epsilon)_+^2/v.                          (7)
```

Indeed `m <= epsilon + E[G 1{G>epsilon}]`. If `G<=B` almost surely, the sharper
bounded-support version is `(m-epsilon)_+/(B-epsilon)`, for `epsilon<B`.
These bound the separated fraction; they do not establish that it is one.

There is also a **sharp marginal-only upper bound on diagonal mass**. Let
`mu_1,mu_2` be the two birth marginal measures and let
`TV(mu_1,mu_2)=sup_A |mu_1(A)-mu_2(A)|`. Assume the marginals admit an ordered
coupling. Then

```text
sup_{ordered couplings} P(S=T) = 1-TV(mu_1,mu_2).             (8)
```

**Proof.** Any diagonal part is a common submeasure of both marginals, so its
mass is at most the overlap `nu=mu_1 wedge mu_2`. Remove this maximal common
submeasure from both marginals. Their cumulative distribution difference is
unchanged and nonnegative, so the equal-mass residuals still admit an ordered
quantile coupling. The two residual measures are mutually singular and hence
put no mass on equality. Combine that residual coupling with the diagonal
coupling of nu. This attains the bound. The zero-residual case is immediate.

For ordinary density marginals the overlap is `integral min(f1,f2)`, but the
measure formulation includes atoms. Also `q(t)=F_1(t)-F_2(t)`, giving the weaker
bound `P(S=T)<=1-sup_t q(t)` directly from rank-one occupancy.

Consequently, overlapping limiting birth marginals generally leave a real
diagonal ambiguity. The new two-time statistic removes information that cannot
be restored by making the same static curves arbitrarily precise. Equation
(8) is an upper bound over all ordered couplings, not a claim that percolation
can realise every such coupling.

Apply this observation to the limiting marginals with care: finite-n total
variation need not converge under weak convergence. In the counterexample,
the finite marginals have disjoint supports and (8) is zero at every n, but
their limit overlap is `1/2`. Thus a finite marginal-overlap bound of zero
does not repair the missing uniform near-diagonal control either.

## 4. Paired permutation archives: exact continuous-gap conversion

For iid continuous Uniform site labels, their random ordering is independent
of their sorted values. If the ranks are determined by that ordering, then
conditional on `J1=i,J2=j`,

```text
(T1, T2-T1, 1-T2) ~ Dirichlet(i, j-i, N+1-j),               (9)
```

with a zero middle spacing when `i=j`. In particular, for `d>=1`,

```text
T2-T1 | D=d ~ Beta(d,N+1-d),
E[T2-T1|d]  = d/(N+1),
Var[T2-T1|d] = d*(N+1-d)/[(N+1)^2*(N+2)].                  (10)
```

For an affine scaled gap `G=a_N*(T2-T1)`, let `h_d=P(D=d)` and
`x=min(delta/a_N,1)`. With `I_x(alpha,beta)` the regularised beta CDF,

```text
P(G<=delta) = h_0 + sum_{d>=1} h_d I_x(d,N+1-d),

Z(delta) = h_0 + sum_{d>=1} h_d [
    I_x(d,N+1-d)
    - (a_N/delta) * d/(N+1) * I_x(d+1,N+1-d)
].                                                        (11)
```

This is exact conditional integration of label spacings, not new samples and
not a Gaussian approximation. A D histogram suffices for (11); a full `(J1,J2)`
histogram additionally retains location and midpoint-gap dependence. A pair of
separate marginal histograms, or a few joint moments, does not supply `h_d`.
Non-uniform labels require their known CDF transformation. A conditional or
biased permutation sampler requires its sampling weights and cannot silently
be read as an iid permutation ensemble.

For a nonlinear canonical clock `b_N(p)`, the gap
`b_N(T1)-b_N(T2)` needs the joint `(i,j)` law and the joint spacing integral in
(9), not only D. An estimated clock introduces shared estimation error; either
propagate it within the same batches or use a separately calibrated clock.

### When insertion and label gaps have the same limit

Write `G_tilde=a_N*D/(N+1)`. Conditional on D,

```text
E[G-G_tilde|D] = 0,
E[(G-G_tilde)^2|D]
    = a_N*G_tilde*(1-D/(N+1))/(N+2).                       (12)
```

If `G_tilde` is tight and `a_N/N -> 0`, then `G-G_tilde -> 0` in probability:
restrict to `G_tilde<=K`, apply conditional Chebyshev to (12), then let K
increase. A bounded mean yields the stronger mean-square bound
`E[(G-G_tilde)^2] <= a_N E[G_tilde]/(N+2)`.

Thus the two gap variables have the same subsequential laws under these
conditions, including their limiting diagonal atom, even though their finite
CDFs differ. Uniform closeness of CDFs at a shrinking delta is **not** implied;
use (11) for a finite-size near-diagonal report.

For an `L x L` lattice (`N=L^2`), the proposed thermal scaling is
`a_N=L^(3/4)`, hence the insertion scale is `(N+1)/a_N`, asymptotic to
`L^(5/4)`. The spacing condition in (12) is then favourable. For centred birth
**locations**, each order statistic has variance at most `1/[4(N+2)]`, so
`a_N=o(sqrt(N))` also makes the two individual clock errors vanish; this
stronger condition holds for the proposed scaling as well.

The exponent `3/4` is a square-site universality hypothesis/control here, not
a theorem supplied by this note. A triangular comparison should use a specified
near-critical normalisation and torus modulus; known triangular ingredients do
not by themselves prove the required torus hitting-time continuity. The exact
identities (1)--(12) do not depend on choosing this exponent.

## 5. A small prospective comparison, not another permission system

Use one paired archive or one independently seeded production block, retain
aligned batch identities, and report the same summaries at every size. The
following are competing interpretations, not claims already established.

| Question | Prespecified readout | Result that changes the mechanism picture |
|---|---|---|
| Separate macroscopic births? | `Z_L(delta)` and beta-mixture `P(G<=delta)` at `delta=1/2,1/4,1/8,1/16`, then smaller fixed resolutions if informative | A size-stable positive small-delta floor challenges simple-birth noncoalescence; decreasing values support it without proving the double limit. |
| Microscopic direct fusion versus aggregation? | `h_0` alongside `Z_L(delta)-h_0` on the same block | Falling `h_0` with persistent soft mass identifies the aggregation alternative that an arm count alone misses. |
| Clock or location escape rather than aggregation? | Birth marginal quantiles and tails, `G_tilde` quantiles; raw thermal scale plus a separately specified canonical clock when available | Escaping locations/gaps invalidate an asserted tight pair limit; clock-dependent near-diagonal conclusions call for clock control, not a new fitted exponent. |
| What has one-time analysis already told us? | `E D=sum q_1(k)` and marginal overlap bound, versus the paired `Z` curve | Equal static summaries but distinct `Z` directly isolate process information, as in the counterexample. |

These rows can be computed together; they are not sequential approval gates.
The numerical delta grid is a reporting convention, not a universal resolution
or significance boundary. Fix the clock and grid before inspecting a new block.
Save the cross-delta/cross-readout batch covariance. Shared random streams or
reanalysed archives remain one dependency group; do not add the several rows as
independent evidence. If the paired archive is missing, report the marginal
bound (8) and the precise missing object instead of inventing a copula.

No GPU or massive trajectory storage is intrinsic to this test: a joint birth
histogram, batch IDs and optional direction marks are sufficient. The next
theory target is a **uniform near-diagonal bound** for this paired object,
potentially through multi-time pivotal geometry. A single-site six-arm bound
only addresses its direct-jump part unless an additional argument controls
aggregation of separate events.

## 6. Executed check and claim boundary

`python3 scripts/check_birth_gap_counterexample.py` performs one small exact
Fraction calculation: equal static marginals in the two counterexample
families; the discrete overlap/curvature identities; and (11) against direct
integration of integer-parameter beta density polynomials. It enumerates no
percolation configurations and draws no Monte Carlo samples.

This note adds the counterexample, occupation-defect criterion, sharp
marginal-overlap bound, and exact finite archive interface. It does not provide
a measured large-size diagonal mass, a new percolation exponent, a proof of
canonical copula universality, or an identification of the original-U operator.
