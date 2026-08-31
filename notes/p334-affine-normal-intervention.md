# A direct independent experiment for the frozen source-normal response

Use a pair of **affine probability laws**, rather than estimating a
derivative by a small exponential tilt. They preserve every existing
joint-safe degree-class mass and the means of both original density scores
at finite amplitude. Their weighted difference equals the archived normal
response exactly, without a Taylor remainder.

## Freeze the source and population before drawing new prefixes

Fix the geometry/size, original checkpoint rule, own source o=f or s,
receiving observer, new independent RNG/counter block, sample budget, and
decision rule. The minimal unit is one predeclared (N,o) with two jointly
recorded outcomes: normalized birth center C and A(p_ref),
`p_ref=.59274605079`. Extra sizes/receivers require predeclared weights or
multiplicity allocation; they must not be selected using the new outcomes.

Generate prefixes from the same original population as the discovery
archive. Keep its full population denominator: original non00 prefixes
have zero normal contribution and need no intervention suffixes. Generating
only00 prefixes and silently removing their population weight changes the
estimand. The old archive supplies the prediction, not additional samples
for the independent test.

At each new00 prefix compute the *already defined*, response-free census:
the two first scores s=(s_f,s_s), own second density score T_oo, and

```
G=E_0[s s^T],       alpha=E_0[T_oo s^T] G^(-1),
phi=T_oo-alpha s.
```

Expectations E_0 are uniform over that prefix's vacant labels. This is the
frozen source-normal rule from `notes/p334-source-normal-curvature.md`,
not another descriptor search or a fit to new responses. It gives

\[
 E_0[\phi\mid A_a]=0,\qquad E_0[\phi s_f]=E_0[\phi s_s]=0,
\]

with phi=0 outside the joint-safe classes. Invertibility in the discovery
prefixes is not an assumption about every independent prefix. Freeze an
exact-rank continuation: if G has rank1, project on the first nonzero score
in the declared f,s order; if rank0, both own second score and phi vanish.
This projects on the same actual source span with no ridge or floating
eigenvalue cutoff. It preserves the class and orthogonality identities.

## Two finite laws, exact conservation, exact target

For B(Z)=max_u |phi(Z,u)|>0 define

\[
 \boxed{q_\pm(u\mid Z)={1\over d}
                       \left[1\pm{\phi(Z,u)\over B(Z)}\right].}
\]

They are nonnegative and normalized. Each class has exactly its original
mass because its phi sum is zero; outside probabilities remain1/d.
Consequently the entire immediate rank/Euler joint distribution of both
geometries is conserved. Furthermore

\[
 E_{q_+}s_i=E_{q_-}s_i=E_0s_i=0,\quad i=f,s.
\]

This finite-score preservation is stronger than the first-order
orthogonality of an exponential normal tilt. It does **not** assert that
each class's score mean, the unweighted raw loop marks, or score covariance
is conserved. Those stronger statements do not follow from the projection.

At B=0 record an exact zero contribution and skip both intervention tails.
Using `max(1,B)` is also mathematically valid, but changes efficiency; the
minimal choice B=max|phi| is the proposed frozen rule. Do not switch the
normalization in response to observed outcomes.

After each selected label use the original uniform suffix law. For every
fixed observer F,

\[
 {B(Z)\over2}\{E_{q_+}F-E_{q_-}F\}
       =E_0[\phi F\mid Z].
\]

Thus `D_F(Z)=B(Z)(F_plus-F_minus)/2`, with any valid paired coupling of
the two branches, directly estimates the same full-population normal
response as the discovery calculation. **B stays inside the population
average.** An unweighted branch difference averages H_F(Z)/B(Z), a
different target whose sign need not match the archived aggregate sign.

There is also a parameter-free mean identity:
`(E_q+ F+E_q- F)/2=E_0 F`, since `(q_++q_-)/2=q_0`. Pooling the two arms
without their contrast therefore restores the unperturbed mean.

## A small paired sampler with exact zero differences on its shared part

Let `rho(Z)=E_0|phi|/B`, the total-variation distance between q+ and q-.
A maximal same-label coupling can be sampled directly:

- With probability1-rho, choose the same label in both arms with mass
  proportional to `1-|phi(u)|/B`.
- With probability rho, choose the plus label proportional to phi_+ and
  the minus label proportional to (-phi)_+.

The equal total positive/negative phi mass proves the required marginals.
Draw one independent uniform permutation of the remaining original labels
and remove the forced first label separately for each arm. Each resulting
suffix is uniform; the two geometries share each arm's complete sequence.
When the forced labels agree, reuse exactly the same suffix and obtain
`D_F=0` samplewise. No variance-optimality claim is needed.

If F has range length R_F, this coupling gives the useful conditional bound
`E[D_F^2|Z] <= B R_F^2 E_0|phi|/4`. The range lengths are1 for normalized
C and2 for A. This is a sampling-cost bound, not a claim about the actual
variance or sufficient production size.

## The frozen competing predictions and a finite stop rule

The strong closure model is, for each prefix and each primary observer,

```
M0: m_F(Z,u)=c_a(Z)+b_f(Z)s_f(Z,u)+b_s(Z)s_s(Z,u),
```

with arbitrary class intercepts and a common pair of slopes across classes.
It predicts **D_C=D_A=0 in expectation at every prefix**, hence zero
population response, exactly at the finite policies above. Class-dependent
slopes are not covered by this null and can produce normal response.

The archived alternative predicts positive own-center and negative own-A
weighted population response on a new source block. Before generation,
freeze whether the test is sign-only or also requires minimum transported
magnitudes `(delta_C,delta_A)>0`. One transparent practical proposal is
one half of the selected archive's point magnitude; that fraction is a
declared replication target, not a mathematical constant. Then

```
M1: H_C >= delta_C,    H_A <= -delta_A.
```

Use one fixed primary analysis at the predeclared maximum independent
prefix count. Construct simultaneous uncertainty intervals for the two
B-weighted responses using the paired prefix/batch vectors. For example,
allocate total alpha=.01 across the two coordinates; if several geometry
units are separately claimed, include them in the allocation. An early
look is allowed only with its timing and error allocation frozen in advance.
No continued sampling until a preferred sign appears.

For the resulting joint region `[L_C,U_C] x [L_A,U_A]`:

| Observed region | Decision and stop |
|---|---|
| L_C>0 and U_A<0 | Replicate the directional prediction and reject M0; the minimum-magnitude version is established only if L_C>=delta_C and U_A<=-delta_A. |
| U_C<delta_C or L_A>-delta_A | Falsify the frozen minimum-magnitude M1; independently report whether the origin is excluded. |
| Origin excluded, but archived directions not established | Reject the tested strong closure without declaring M1 or a new field count. |
| Remaining overlap at the maximum budget | Stop as unresolved at the stated precision; do not retune phi, B, outcomes or thresholds. |

For a sign-only alternative there is no positive effect-size floor; a
near-zero result cannot by itself falsify an arbitrarily small signed
effect. Conversely, a zero population response cannot prove the prefixwise
closure because normal responses of different prefixes can cancel. M0
and the specified M1 are disjoint predictions, not an exhaustive list of
all mechanisms. A result can reject both.

## What a successful independent intervention would mean

It would show an operational finite-label response outside the two original
score span under a rule that leaves both immediate rank/Euler laws and
both score means unchanged. It would independently test the archived
normal C-positive/A-negative prediction. Those two outcomes share birth
clocks and paired sampling; they are not two independent replications.

This does not establish that the source causes the unperturbed global
Matching-One anomaly, identify a continuum operator, count fields, or
invalidate a successful first-Jacobian predictor `J=B G`. The old source
law is exactly the mixture of the two new arms. The experiment measures
the response to a deliberately selected probability change, not the
contribution of that direction to the untouched mean.

Theory receipt: follows discovery2c3a5ca2 and its f97c227c independent-test
handoff. No new prefix, suffix, descriptor, model fit, numerical test or
server task was run here. The root owns the frozen production protocol,
generation, covariance analysis and publication of its actual outcome.
