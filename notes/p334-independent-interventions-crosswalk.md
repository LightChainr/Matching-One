# Two independent P334 interventions constrain different pieces of the label mean

**The results coexist exactly.** The frozen four-contact experiment rejects
both near-zero and near-old-value predictions for one *across-prefix loading
of first-score responses*. The source-normal experiment detects a different
piece of the conditional label mean, orthogonal to those first scores. Its
positive result rejects their promotion to a complete label-mean description.
The observed contact ratio near one half is not a new half-amplitude model.

## Exact crosswalk at each complete prefix

Use uniform next-label expectation conditional on prefix Z. Let
`m_F(Z,u)=E_suffix[F|Z,u]`, and H=(H_f,H_s) be the two frozen first density
scores. Each H has zero mean in every joint-safe degree class a and is zero
on unsafe labels. Define, on safe labels,

```text
m_F^c = m_F − E[m_F|a,Z],
G = E[H H^T|Z],           tau_F = E[H m_F|Z],
b_F = G^+ tau_F,          r_F = m_F^c − H^T b_F.
```

Set m_F^c=r_F=0 outside the safe set; outside m_F need not be zero.
G^+ means the inverse on the exact Gram range, including rank1/0. If
v is in ker G, then v^T H=0 almost surely and v^T tau_F=0. Hence tau_F is
in range G, proving without an invertibility assumption

```text
tau_F = G b_F,    E[H r_F|Z]=0,    E[r_F|a,Z]=0.
```

For the frozen own-axis normal score phi_oo, which is class-centered,
zero outside the safe set, and orthogonal to both H components,

```text
nu_oF(Z) = E[phi_oo m_F|Z] = E[phi_oo r_F|Z].
```

Thus the first response sees the projection coefficients b_F through G;
the normal response sees only r_F. Changing r_F within the orthogonal
complement cannot change tau_F. Changing b_F cannot change nu_oF. This is
an exact decomposition, not a claim that b and r are independent random
variables across prefixes, or that the identity tau=Gb is itself a
predictive model. The coefficient b_F is also **not** the four-contact
regression coefficient frozen from the training archive.

For receiver o, that fixed predictor estimates the own-source component
`tau_oF(Z)=(G b_F)_o` from four prefix features. Its residual functional is

```text
R_o = pi00 * [2 Cov_00(mu_C, (G b_C)_o − rhat_oC)
             − 0.5 Cov_00(mu_W, (G b_W)_o − rhat_oW)].
```

Here mu_F=E[F|Z], C=(K1+K2)/[2(N+1)], and W=(K2−K1)/(N+1).
R is a signed cross-prefix covariance projection; it is neither a norm of
r_F nor the mean E[nu_oF]. There is no identity converting an R ratio into
a normal-response ratio. Since H and r_F both have zero label mean, mu_F
depends on class means and the unsafe-label mean, not on a hidden average
of r_F. In particular, failure of a forecast for R and
success of a forecast for the mean normal response impose no conflicting
conditions on the label mean.

## What the two completed independent experiments found

| Experiment | Perturbation and primary functional | Fresh population and result |
|---|---|---|
| Frozen four-contact forecast | `q±=(1±H_o/8)/d`; exact unit response tau_oF, then paired-receiver R above | 300,000 prefixes/N; 60 new batches/N. R_new/R_old=.4988857 at N325 and .5169035 at N425. |
| Normal label-mean closure | `q±=(1±phi_oo/B(Z))/d`; prefix-weighted `B(F+−F−)/2=nu_oF` in expectation, then equal four-own-C mean T | 500,000 prefixes/N; 20 new batches/N. T=(3.0852005663±.3918738407)×10^-8, one new-batch SE. |

Both retain non00 prefixes as zero contributions with the full-population
denominator, but their validation streams are separate new RNG domains and
their policies and functionals differ. Within each experiment the receiving
geometries and derived clock/readout coordinates retain shared-batch
covariance. Both use frozen discovery information, not pooled old/new
validation data.

For the contact experiment the intervals from the prespecified rule are
`[.4360616,.5617098]` and `[.4506760,.5831311]`. Both exclude the C0 band
`[−.25,.25]` and C1 band `[.75,1.25]`. This rejects the fixed claims that this
loading vanishes or transports within25% of its frozen old point value.
It does not turn C0 into a claim of complete four-feature response closure.
Source: completed
[`14b2c98e:REPORT.md`](https://github.com/LightChainr/Matching-One/blob/14b2c98ed3a252a2fe79ce5e124d9484b23a264f/experiments/p334-prospective-intervention-20260831/REPORT.md),
with its frozen F0 contract4b3c21b7 and separate60-batch inference.

For the normal experiment, T±3SE is `[1.9095790443,4.2608220883]×10^-8`.
Its lower endpoint exceeds the predeclared practical-zero threshold10^-8;
the frozen archive forecast3.6565×10^-8 lies inside that interval without
amplitude refitting. A complete safe-label mean
`m_C=c_a(Z)+b_f(Z)H_f+b_s(Z)H_s`, with a common slope pair across classes,
would have r_C=0 and all normal responses zero. This is the stronger model
rejected by the planned rule. A relation predicting only the first Jacobian
does not impose r_C=0. Source:
[`5d19fe79:normal-intervention result`](https://github.com/LightChainr/Matching-One/blob/5d19fe79497050c59b5e43365c0534b4cd333f73/notes/p334-independent-normal-intervention-result.md),
score1164ba91, frozen decision43079652/bc0a18c2.

## Unified physical meaning and one boundary

The finite interventions establish future birth sensitivity while preserving
the specified instantaneous rank/Euler class law. The first experiment limits
transfer of a frozen predictor for the source-visible coefficient field;
the second shows a nonzero source-orthogonal component of the future center
mean. Class-dependent susceptibilities or higher within-class structure can
produce the latter without contradicting a successful first-response law.

These are two distinct completed tests, not votes for one field or a common
effect size. The contact ratio is conditional on frozen training parameters;
it does not isolate training error from model form, and receives no post-hoc
half-amplitude rescue. The normal3SE rule is the frozen fixed-batch diagnostic,
not exact finite-sample coverage or confirmation of its forecast amplitude.
Neither result identifies the source of the unperturbed global anomaly, counts
fields, or selects a complete microscopic mechanism. No extra descriptor,
fit, sampling, raw replay, combined significance test or production proposal
was added in this crosswalk.
