# The threshold law's finite-size change is not center and width, and ~99% of what is left is one direction

**Date:** 2026-09-06
**Ticket:** #582 (first data pass)
**Claim level:** C2 — reanalysis of committed productions, no new samples
**Artifact:** `results/wasserstein-shape-flow/latest.json`
**Scripts:** `scripts/threshold_quantile_lineage.py`, `scripts/score_wasserstein_shape_flow.py`

#582 asks whether the finite-size motion of the whole threshold law is only a
center/width drift or a genuine change of shape, and — if the latter — whether
the shape flow is one-dimensional, two-dimensional, or broad and rotating. This
is the first pass on real data. Five transitions across three lineages and two
independent production campaigns, no new sampling.

## The pipeline

`F_N(p)` is the committed channel of `analysis/threshold_histogram_profile_contract.json`
— a transition at the k-th occupied site is `Beta(k, N+1−k)`, and the mixture is
the equal-weight average of the `K_minus` and `K_plus` rank laws. The quantile
function is inverted per **batch**, so an aligned delete-one that removes batch
`b` from every level and both orientations at once gives the covariance of the
whole grid. The affine tangent `span{1, Q_N}` is a two-column basis for
`projective_inference.subspace_residual` — #582's shape flow is #579's statistic
with `dim(V) = 2`. Nine frozen levels leave 7 df.

Production does not use the exact rational path: at N=325 the degree-N polynomial
carries `N·C(N−1,k−1)` and Sturm isolation is unaffordable. The identity

```text
F(p) = Σ_j B(N, j, p) G(j),      G(j) = P(K ≤ j)
```

collapses the mixture into one binomial pass, anchored at the mode so the
`(1−p)^N` underflow that killed the N=1300 pilot cannot recur by construction.
Checked against the exact rational path at N=4, 12 and 20: CDF agrees to
`4e−16`, and every bisection quantile lands **inside** the exact Sturm bracket.

## The correction I had to make to my own first pass

The first run averaged the two orientations with equal weight. That is wrong, and
wrong in a way that nearly produced a finding.

Two Gaussian integers of the same norm do not have opposite `cos 4θ`. The
equal-weight net residue is

```text
N=  65 +0.1972     N=  85 -0.1560     N= 145 -0.0137
N= 130 -0.1972     N= 170 +0.1560     N= 290 +0.0137
N= 325 +0.5937     N= 425 -0.4697
```

It **alternates in sign along each lineage** and has **opposite sign between the
two Gaussian lineages** — exactly the pattern of the lineage split the first run
showed in the shape residual. So the law is now reconstructed twice: the naive
average, and the combination `w₁Q₁ + w₂Q₂` with `w₁c₁ + w₂c₂ = 0`, which removes
spin 4 exactly to first order.

At N=325 and N=425 both orientations carry the same sign of `cos 4θ`, so those
weights leave `[0,1]` and the correction is an extrapolation rather than a
mixture. That is recorded per size rather than hidden.

**The correction is real and does not explain the split.** It moves the median by
0.1σ to 5.0σ, and the angle between the two lineages' shape directions goes from
47.0° to 42.4°. So the worry was legitimate and the finding survives it.

## Control 1 — the pipeline does not manufacture shape flow

Two **disjoint halves of the same production**, scored as if they were a
transition. Nothing changes between them, so the shape flow must be zero.

| N | D | df | p | ‖v‖ |
|---:|---:|---:|---:|---:|
| 65 | 4.93 | 7 | 0.669 | 5.8e−05 |
| 145 | 8.55 | 7 | 0.287 | 9.3e−06 |
| 290 | 4.13 | 7 | 0.765 | 1.9e−05 |
| 325 | 2.08 | 7 | 0.955 | 6.7e−06 |

Consistent with zero at every size, and the half-to-half displacement is a
thousandth of a real transition's. This is the control that makes a `χ²` of `10⁵`
below mean something.

## Result 1 — W0 is dead

| transition | m | D | df | ‖residual‖/‖v‖ |
|---|---:|---:|---:|---:|
| 65→130 | 2.0 | 435,744 | 7 | 2.4% |
| 130→325 | 2.5 | 192,537 | 7 | 1.4% |
| 85→170 | 2.0 | 292,805 | 7 | 3.2% |
| 170→425 | 2.5 | 154,877 | 7 | 2.2% |
| 145→290 | 2.0 | 112,067 | 7 | 1.8% |

The threshold law's finite-size change is **not** a center shift plus a width
change. The shape part is small — one to three percent of the displacement — and
overwhelmingly resolved, because 100M–500M samples put the quantile standard
errors at `3e−06` against a residual of `7e−04`.

Two sanity anchors: `Q(0.5) = 0.59274` at every one of the eight sizes, to a few
parts in `10⁶` — the intrinsic quantile center transfers, which is #101/#236's
result recovered here for free — and the width `Q(0.9) − Q(0.1)` falls smoothly
from 0.243 at N=65 to 0.123 at N=425.

## Result 2 — one frozen direction carries ~99% of it, and it transfers

A shape generator extracted from **four** transitions and never refitted on the
fifth, with only its amplitude free there:

| held out | D affine | rank-1 removes | rank-2 removes | random median | best of 200 random |
|---|---:|---:|---:|---:|---:|
| 65→130 | 435,744 | **99.08%** | 99.75% | 1.14% | 19.48% |
| 130→325 | 192,537 | **99.77%** | 99.92% | 0.79% | 17.26% |
| 85→170 | 292,805 | **99.82%** | 99.98% | 3.34% | 25.35% |
| 170→425 | 154,877 | **99.29%** | 99.66% | 5.88% | 24.79% |
| 145→290 | 112,067 | **99.76%** | 99.99% | 3.72% | 22.22% |

The random-direction control is what makes this a result. A direction that was
never fitted removes 1–6% of the `χ²`; the *best* of 200 random draws removes
17–28%; the frozen trained direction removes 99.1–99.8%. On 145→290 the generator
was trained entirely on the other two lineages, different productions and
different seeds, and still removes 99.76%.

Grid robustness — the frozen grid was fixed before any lineage was loaded, and
rerunning everything on three others changes nothing:

```text
frozen 9 deciles   df=7   rank-1 removed 99.08%-99.82%
narrow 0.25-0.75   df=6   rank-1 removed 99.09%-99.81%
wide  0.05-0.95    df=7   rank-1 removed 99.12%-99.80%
fine  15 levels    df=9   rank-1 removed 99.04%-99.84%
```

## Result 3 — and it still does not close, which is the interesting part

After the frozen rank-1 generator the residual is still enormously significant,
and rank 2 does not fix that either. So on #582's decision table this is **not**
W0, **not** W1 and **not** W2 — and just as clearly **not** Wbroad. It is a fifth
outcome the table does not have:

> a dominant transferable direction carrying ~99% of the flow, plus a remainder
> that is small, highly resolved, and structured.

The structure in the remainder is visible in the pairwise angles between the
transitions' unit shape residuals:

```text
85->170   vs 170->425     3.9 deg
65->130   vs 145->290    14.1 deg
130->325  vs 85->170     16.6 deg
65->130   vs 85->170     42.4 deg
85->170   vs 145->290    55.7 deg
```

`{85→170, 170→425}` agree with each other to 3.9°; `{65→130, 145→290}` agree to
14.1°; the two groups sit 42–56° apart. This is **not** the spin-4 residue (the
correction barely moves it), **not** a production artifact (85 and 170 come from
the same runs as 65 and 130), and **not** the scale multiplier (each group mixes
m=2 and m=2.5, or two m=2 from different campaigns).

What indexes the second direction is **not determined by five transitions**. It
is not monotone in `N` — the clusters are `{65,145}` and `{85,170}` with 130 in
between — so it is not a smooth finite-size correction in the obvious sense.

## What this means for the decision tree

The directive that opened this work makes #582 the gate: a stable low-dimensional
shape licenses talking about a compact RG state; broad/rotating would mean the
project's low-rank results are observer bandwidth.

The answer is **closer to the first than to the second, and is not either**. There
is a real, transferable, cross-lineage shape generator — that is a genuine
compact object and it survived a held-out test against a random-direction control
by a factor of ~20 over the best of 200 draws. But it does not close the flow,
and the part it misses has reproducible structure rather than noise.

So the honest input to the next step is: **rank ≈ 1 with a resolved correction**,
not "rank 1" and not "broad". #581 can be asked its question — is the dominant
direction bulk-like or ambient-topology-like — because there *is* a dominant
direction to ask about. What #582 cannot yet hand over is a clean rank, and the
place that is decided is the second direction, which needs more base sizes or a
readout that separates them.

## Not established

- **Any mechanism.** A shape direction is a direction. Naming it is #581's job.
- **Freedom from spin 8.** The correction removes spin 4 to first order only, and
  N=325 and N=425 need an extrapolation rather than a mixture to get even that.
- **An exponent for the shape amplitude.** The rank-1 amplitudes on the held-out
  transitions are −2.14, −1.40, −0.95, −0.90, −0.65 (×10⁻³), and three base sizes
  do not determine a power. This analysis deliberately does not fit one.
- **That the cross-size covariance is zero.** It is measured batch by batch: the
  largest per-level correlation is 0.10–0.19 against a noise floor of 0.10, so it
  is dropped on evidence rather than by assumption.
- **Anything about the orientation-difference channel**, which carries the spin-4
  amplitude and is a different observable from the law reconstructed here.
