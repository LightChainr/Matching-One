# A hole reweights the intact topology: the jump-only mixed-U model fails

**New exact decision:** the baseline-reweighting contribution is
`U_st^rw = −4.550327123236791`. Its reduced rational enclosure strictly
excludes zero. Therefore the fixed **weighted-rank-jump-only** model fails
for the original global U on this N50 pair. No source, coefficient, chart
or defect class was selected from this outcome.

The source remains `S*=Ctot+F4+Bvac`, with bulk exp(tS*), and the chart is
`pA=s+(1−s)p, pB=p`. Increasing s means moving toward saturation;
for epsilon=1−s the same contribution is `U_t,epsilon^rw=+4.550327123236791`.

## The fixed operator split

Let O0,S0 denote the intact parent, O1,S1 its origin-A-vacancy version, and
DeltaS=S1−S0. Both are evaluated on the **same free-B configuration**.
With normalized intact expectations and w=exp(tDeltaS),

```
partial_epsilon <O> = 25(1−p) [<w(O1−O0)> + Cov(w,O0)].
```

The second term is identically zero at t=0 for every p. Its t derivative
there is Cov_p(DeltaS,O0), without an extra intact-source covariance:
derivatives of Cov(1,O0) cancel. Thus its mixed mean jet in the s chart is

```
C_O(p)=Cov_p(DeltaS,O0),
g_O(p)=(partial_st <O>)_rw=−25(1−p) C_O(p),
g'_O=25 C_O−25(1−p) C'_O,            O=q,E.
```

The factor25 is translation over all A holes; the explicit derivative of
1−p is retained. Per-geometry normalization is retained through

```
C_O = <S1 O0>−<S0 O0>−(<S1>−<S0>)<O0>.
```

Only `<S1*q0>` and `<S1*E0>` were absent from the committed statistics.
Neither a product of marginal means nor matching two independent
enumerations would supply them. The new paired observer saves precisely
these two integer cross moments, plus configuration counts, at each k.

## Why the mixed contribution is a linear root-complete readout

Put Q=mean(q0), Y=P4(E0), D=Q_p and A=50^(13/8)/2 at the same saturated
pooled root. Because the reweighting s insertion and **all its p
derivatives** vanish at t=0, differentiating the U response operator in t
has no additional product with that zero insertion. The surviving term is
the ordinary complete linear response operator applied to g:

```
U_st^rw/A = gY_p/D − Y_pp*gQ/D² − Y_p*gQ_p/D²
            + Y_p*Q_pp*gQ/D³.
```

In particular the root contribution is `p_st^rw=−gQ/D`; root and slope
motion have not been suppressed. This linearity defines an additive
operator contribution, not a separate empirical population or a fitted gain.

## Exact numerical readout and two transparent decompositions

| Fixed component of U_st | Value |
|---|---:|
| **Baseline reweighting: new primary** | **−4.550327123236791** |
| Weighted rank jump: old total minus new primary | +15.306045530800864 |
| Complete U_st: imported earlier enclosure, not recomputed | +10.755718407564073 |

Reweighting opposes the weighted-jump contribution. The complete Xi and
gain residual R remain the previous results; neither was recalculated or
counted as a new outcome here.

Within the new primary, a nonoverlapping dose/root split gives

| Term | U_st contribution |
|---|---:|
| Thermal p derivative of the covariance | −3.8340753980328643 |
| Explicit derivative of the physical hole dose | −0.7342519482682530 |
| Root relocation, including slope relocation | +0.0180002230643270 |

An **alternative**, not additional, linear split has uncentered cross-moment
difference `+0.8646659704866687` and covariance-centering term
`−5.414993093723459`. Dropping normalization would therefore reverse this
contribution's sign. Full rational term bounds are retained in
[`score.json`](../results/p337-defect-reweighting/score.json).

## What has—and has not—been localized

The hole weight is fixed by `DeltaS=3−2*k_null−ell`. Its reweighting term
can act when rank does not change, but it also contains configurations
where rank changes. This calculation imposes **no ell=0 population gate**.
Consequently it does not measure a rank-preserving share, an event fraction,
or a percentage of a field. What is excluded is exactly the model which
omits the normalized reweighting term after projection to this global U.

This remains the fixed axis/tilted N25-child pair with different Smith
classes, observed through its two N50 parents. No continuum field,
asymptotic exponent, replacement source or new production block follows.

## Frozen input and execution

- Contract: `e6a900d9d644b26278f01c17bdfb6f27f3903b75`.
- Producer/scorer: `db348346`, before any new enumeration or scoring.
- Input repository: `359bde9b`; specified question: overview `7132f0c2`.
- Each parent, (5,5) and (1,7), made one complete paired 2^25 traversal,
  preserving the original quotient/rollback backend. No separate saturated
  pass ran. The new counts took1.79/1.83seconds; compilation, both traversals
  and the single rational scoring pass took3.70seconds overall.
- The saved child-root enclosure was complemented to
  `p0≈0.40733446067177326`; there was no new root search.
- Old intact and defect marginals were reused. The old complete U_st/A
  enclosure was imported only to form the complementary linear term.

The result folder records all input/source hashes and the rational primary
enclosure. These are deterministic computational bounds conditional on the
stated exact graph counts, not confidence intervals. No Monte Carlo, cloud
job, extra fit or test suite was run. This branch is not pushed by this task.
