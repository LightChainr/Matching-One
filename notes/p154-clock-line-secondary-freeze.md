# P154 secondary clock-line decision: fixed before reading fresh outcomes

This fixes the operational comparison of the **two already proposed**
clock fingerprints in `c2828e3430fe1ac7e02fbe0e5ddc0e6a24c99847`.
It does not change the primary experiment, introduce a fourth numerical
rectangle, fit a new source, or collect data. The executor has read the
published production code/contract and verified a live process, but has
not read any fresh P154 source response or prospective score. P334 is a
different completed intervention and is not an input to this comparison.

## Fixed population, source, and hypotheses

- Original experiment freeze: `0820b8d203e2dc534bb883d6fdb4d6d1e0acb11f`.
- Exactly the completed N85 5,000,000 and N340 160,000,000 samples, with
  200 batches per N and all nine planned shards. No partial scoring,
  additional samples, size substitutions, or old/fresh pooling.
- The original lag=1 rank-centered bulk CB+CW source, original geometries,
  canonical thermal differentiation, pooled root and normalized U.
- M10: a01=a11=a01'=a11'=a10'=0 at each specified N's population root;
  a10 may be any real number. M11: a10=a01=a10'=a01'=a11'=0, a11 free.
  The common a00 and a00' may be arbitrary in both. Both hypotheses include
  zero and neither exhausts the physical mechanism space. No amplitude is
  estimated to make the test pass.

The symbols and exact four-clock transmission map are defined in
[the pre-existing theory](p154-birth-clock-transmission-map.md). A failure
rejects these **pure, locally flat first-jet restrictions**; it does not
reject every p-dependent/mixed clock or identify a continuum field.

## Fixed readout and nuisance evaluation

Let C_e(N,m), C_c(N,m) be the entry/completion value gains of mode m,
computed from the same unmarked population baseline jets by the previously
published formula. For each N=85,340 and m=10,11 define

```text
d_(N,m) = [C_e*v_completion − C_c*v_entry] / hypot(C_e,C_c).
```

This signed perpendicular residual is in original-U response units. It is
zero for the specified clock line. We never divide by an observed entry
response or select its sign. The four-coordinate ordering is
`N85.M10, N85.M11, N340.M10, N340.M11`.

Use the official final `PROSPECTIVE_RESULT.json` central entry/completion
values and its already saved 200 delete-one vectors. Do not rerun or
overwrite the official primary score. From the corresponding complete
shards, read **only unmarked sum_q,sum_e** and identities/counts to recover
the baseline root and first/second p-jets. Recompute these baseline gains
inside the same delete-one batch as the official response vector; do not
freeze the noisy old nominal gain table or add its variance as if it were
independent. This estimates the hypothesis's baseline-dependent coefficients,
not a fitted source clock or angular gain.

Use the production contract's bracket [0.55,0.65] and numerical integration
convention. Nonfinite jets, failed bracket, zero gain norm, incomplete or
mismatched batches produce an unscorable result; no changed bracket or
alternative mode is allowed. Retain the N-wise independent sampling domains
and within-N shared geometry/source/baseline dependence.

## Fixed inference and stopping rule

Form paired delete-one factors `sqrt(199/200)*(d_loo−mean(d_loo))` in each
N. Save their full four-coordinate covariance, and a joint factor/covariance
with the official six entry/completion/net coordinates. These are correlated
interpretations of the **same** fresh experiment, never additional evidence
votes. Numerical gates and model coefficients are not revised after reveal.

The secondary family is these four residuals only. Its nominal simultaneous
95% intervals are `d ± z*SE`, with
`z=Phi^−1(1−0.05/(2*4))`. This is Bonferroni with asymptotic normal marginal
intervals, not an exact finite-sample or anytime guarantee. A mode's
two-size restriction is contradicted if either size's interval excludes
zero. Otherwise report **not excluded**, not accepted or a forced winner.
No equivalence margin is introduced and no event-current, Xi, source-alpha
fit, further mode, size, or derivative scan is added to this family.

The official primary decisions remain exactly as delivered. If its three
numerical rectangles fail, they stay failed regardless of this secondary
interpretation. If these two theory restrictions also fail, report both
failures and stop this scorer; do not fit a mixture or a replacement line
to the same block. A secondary null alone never identifies a mechanism.

## Execution receipt

The implementation must be committed before it reads fresh outcomes. Its
receipt records this definition commit, its implementation commit, the
official result commit/hash, the production freeze, all nine input receipts,
and the retained baseline/source dependency groups. Output is a new path;
neither another team's files nor the official score are edited. This is a
single deterministic scoring pass after complete delivery, not permission
to start, restart, extend, or interrupt production.
