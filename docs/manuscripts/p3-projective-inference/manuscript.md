# Test the ray, not the ratio: denominator-free model comparison for finite-size amplitude ladders

**Manuscript draft — P3 of the publication portfolio. Issue [#579](https://github.com/LightChainr/Matching-One/issues/579).**

**Target venue:** *Physical Review E* (Statistical Physics / computational methods). The result is a statistical
method with a percolation application; PRE is where the amplitude-ladder literature this corrects actually lives.

**Status:** complete draft of sections 1–8. Every number below is generated from committed artifacts by
`scripts/p3_manuscript_evidence_table.py` and rendered into [`tables.md`](tables.md); none is hand-typed.
Literature items not yet read in the primary are marked **[LIT]**. See [`README.md`](README.md) for the
section-by-section readiness ledger.

---

## Abstract (draft)

A finite-size scaling model rarely predicts an amplitude; it predicts how an amplitude changes across a ladder of
geometries. Comparing such a model with data is therefore almost always done by forming a ratio of two measured
amplitudes and testing it against the predicted ratio. We argue that this step is a mistake of type, not of
technique, and that it has a consequence which is easy to miss: the practitioner who nominates a denominator
inherits that rung's conditioning, and the natural remedy — dropping a rung you distrust — removes exactly the
information that would have told you whether the distrust was calibrated.

A model that predicts proportions predicts a **ray** in the response space, and a family of such models predicts a
**subspace**. The natural statistic is then the covariance-weighted distance from the measured response vector to
that ray or subspace, `D = min_a (y − Va)ᵀ S⁺ (y − Va)`, which is asymptotically `χ²` on `rank(S) − dim(V)`
degrees of freedom and nominates no coordinate as a denominator. We prove and verify numerically that for a
two-entry response and a one-dimensional model this statistic is *exactly* Fieller's `z` squared, so the
construction is a strict generalization of Fieller's theorem to arbitrary response dimension and arbitrary model
dimension rather than a competitor to it. Within the same framework, model *classes* are annihilated by explicit
linear functionals of the response: on rungs `r = 1, 2, 4` the second divided difference `f[1,2,4]` returns
exactly 1 on `r²` and exactly 0 on any line, with no denominator and no matrix inverse anywhere in it.

We apply this to a prospectively frozen three-rung aspect-ratio ladder for the spin-4 amplitude of a
matching-odd percolation readout at 580 sites, against eight competing modular and non-modular laws. Restricted
to the two rungs the frozen design actually used, the projective statistic reproduces Fieller's `z²` to a largest
relative deviation of `1.5·10⁻¹⁵`, so every subsequent change is attributable to the third rung and not to a
change of statistic. Using all three rungs, two verdicts flip from *compatible* to *excluded* at 7.1σ and 7.0σ,
one exclusion weakens, and seven of eight competitors are excluded at 7.0–11.3σ.

The surviving structure is not a winner. The measured second divided difference is negative, `−4.66·10⁻⁴ ±
1.53·10⁻⁴` (`z = −3.05`), while every competitor predicts it to be exactly zero or strictly positive: the
response is concave in `r` and no candidate law is. Reconciling the dropped middle rung with any competitor
requires a spin-8 to spin-4 amplitude ratio between 7.7 and 785, against the `≪ 1` the frozen design assumed
when it dropped that rung — and the one competitor needing a small ratio, 0.4, is excluded at 9.2σ by the two
clean rungs alone. The design's justification for discarding a rung and the data on that rung are therefore in
quantitative contradiction, a contradiction the ratio test could not express because it lives in the direction
the ratio quotients away.

We conclude with the design consequence: the amplitude the analysis assumed must be measured, which requires
three lattice orientations per family rather than two, and we give the smallest geometry where that is possible.

---

## 1. Introduction

### 1.1 Models of finite-size scaling predict proportions

Very few finite-size scaling predictions are predictions about a number. Overall amplitudes carry
non-universal lattice normalization, and a theory that fixes them is rare. What a theory usually fixes is how an
amplitude *moves* when the geometry changes: between aspect ratios, between boundary conditions, between
lattices in a family. So the prediction that reaches the data is a proportion, and the comparison is made
against a measured proportion.

The received procedure is direct. Measure the amplitude at two geometries, `X` and `Y`; form `R̂ = Y/X`;
propagate the errors; test `(R̂ − R₀)/SE(R̂)` against a normal. The first two steps are unimpeachable. The third
is where a real problem enters, and the fourth is where it becomes invisible.

### 1.2 Two failures, of different kinds

The first failure is distributional and well known. The ratio of two jointly normal quantities is not normal;
when the denominator has non-negligible mass near zero the ratio is heavy-tailed, and in the limiting case a
ratio of independent centred normals is Cauchy, which has no mean. Gleser and Hwang [3] proved the sharp form:
in errors-in-variables and ratio problems, no confidence set for the ratio parameter with finite expected
diameter can have correct coverage. Any procedure that reports a finite symmetric interval for a ratio is
therefore either exploiting a strong denominator or misreporting its coverage.

Fieller [1] solved this in 1954, and the solution is exactly the right one. Instead of studentizing the ratio,
test the linear contrast

```text
Y − R₀ X,
```

whose distribution is normal under the central limit theorem that applies to `X` and `Y` themselves, and invert
over `R₀`. The resulting confidence set is a possibly unbounded region — an interval, the complement of an
interval, or the whole line — and its unboundedness is not a defect but the honest statement Gleser and Hwang
require. Fieller's theorem is not obscure, and part of the point of this paper is that adopting it is necessary
and not sufficient.

The second failure is one of *conditioning*, and Fieller does not touch it. Nominating `X` as a denominator
makes every statement about the model conditional on how well `X` was resolved. Two geometries measured to 1%
each will produce a ratio known to about 1.4% if both are far from zero, and a ratio known to nothing at all if
one of them is not — and *which* one is a choice the analyst made, not a fact about the physics. Worse, the
choice is usually made for a reason that has nothing to do with conditioning: the smallest geometry is cheapest,
so it is measured first and becomes the denominator by default, and the smallest geometry is also frequently
the one where the amplitude is smallest.

### 1.3 The design pathology

The two failures compound through a specific and, we suspect, common design move.

Suppose the ladder has three rungs and the analyst has a reason to distrust one of them — a known systematic
that contaminates that rung and not the others. The natural response, inside the ratio framework, is to keep
that rung out of the decision: form the ratio from the two clean rungs, report the third for information. This
is scrupulous. It is also self-sealing. The size of the systematic was never measured; it was bounded by a
plausibility argument. And the rung that was dropped is the only place where the data could have contradicted
that bound.

We will show that this is exactly what happened in the analysis we reanalyse here, and that the contradiction is
large: the dropped rung requires the systematic to be between 8 and 785 times its assumed bound.

### 1.4 The reframing

The fix is not a better ratio estimator. It is to notice that the ratio was never the object.

A model that predicts proportions across `n` geometries predicts a **ray** through the origin of `Rⁿ` — a
direction, with the overall amplitude unspecified. A one-parameter family of such models predicts a
**subspace**. Neither object requires a coordinate to be singled out. The measured response is a vector `y ∈ Rⁿ`
with an estimated covariance `S`, and the question "is the model compatible with the data" has an answer that
does not mention any ratio:

```text
D = min_a (y − V a)ᵀ S⁺ (y − V a),
```

the covariance-weighted squared distance from `y` to the model's span, which is asymptotically `χ²` on
`rank(S) − dim(V)` degrees of freedom. This is the statistic we advocate. It is not new mathematics — it is the
generalized-inverse quadratic form of Rao and Mitra [5] applied to a linear model without intercept — and its
two-dimensional, one-ray case is Fieller's theorem. What is new here is the observation that this is the *right
default* for finite-size amplitude comparison, and a worked demonstration that the default currently in use
inverts a scientific conclusion.

### 1.5 What this paper claims

1. **A generalization, not an alternative.** The projective statistic reduces to Fieller's `z²` exactly on two
   entries and one ray (§2.3), verified analytically and numerically to `10⁻¹⁴` across six orders of magnitude in
   the predicted ratio, and to `1.5·10⁻¹⁵` on the real data of §4.
2. **Model classes have annihilators.** Whole families of competitors are killed by a single linear functional of
   the response, with no denominator and no matrix inverse, and therefore with no conditioning pathology at all
   (§2.5).
3. **The difference is not cosmetic.** On a prospectively frozen ladder, the two procedures return different
   scientific verdicts for two of eight competitors, and the discrepancy is entirely attributable to a rung the
   ratio framework had a reason to drop (§4).
4. **The residual is informative even when nothing survives.** The response's curvature has a sign no competitor
   can produce, which converts "none of our models fit" into a specific, costed measurement (§6).

We are explicit about what this is not. It is not a claim that the physical law is any of the eight competitors,
nor that it is none of them; §5 sets out what the reanalysis leaves open, including one verdict that is not
stable against a covariance entry the original artifact did not store.

---

## 2. The statistic

### 2.1 A ray, and a subspace

Let `y ∈ Rⁿ` be a measured response across `n` geometries — `n` amplitudes, one per rung of a ladder — and let
`S` be a consistent estimate of its covariance. Batch means, block bootstrap and delete-one jackknife all
produce such a pair; what matters is that `S` estimates the covariance of the *same* random object across all
`n` coordinates, which in practice requires the rungs to share a random stream. We return to this in §7.

A model that predicts the response up to an overall constant supplies a vector `v ∈ Rⁿ` and asserts

```text
E[y] = a v,   a ∈ R unknown.
```

Its prediction is the ray `{a v}`, equivalently the line `span(v)`. A model with `k` free amplitudes supplies a
matrix `V ∈ R^{n×k}` and predicts `E[y] ∈ span(V)`.

Nothing in this statement distinguishes a coordinate. That is the entire content of the reframing: the
prediction is about a *direction*, and directions do not have denominators.

### 2.2 The residual distance

Define

```text
D(V) = min_{a ∈ R^k} (y − V a)ᵀ S⁺ (y − V a),
```

where `S⁺` is the Moore–Penrose pseudo-inverse. Under `E[y] ∈ span(V)` and asymptotic normality of `y`,

```text
D(V) ~ χ²_{ν},   ν = rank(S) − dim(span(V)),
```

by the standard result on quadratic forms in normal variables with a generalized inverse [5]. The pseudo-inverse
rather than the inverse is not fastidiousness: `S` is routinely rank-deficient in exactly the applications this
targets, because the number of rungs can exceed the number of independent random blocks, and because exact
linear constraints among rungs (shared normalizations, common seeds) can make `S` singular by construction. The
statistic remains well defined and the degrees of freedom adjust.

Our implementation (`scripts/projective_inference.py`) forms `S⁺` from a symmetric eigendecomposition in
50-digit arithmetic, discarding eigenvalues below `10⁻¹²` times the largest, and reports the discarded count, the
retained rank and the condition number alongside `D`. Reporting the rank is not optional. A `D` computed at a
rank the analyst did not expect is a different test from the one they intended, and silently inverting a nearly
singular `S` is the same class of error as silently dividing by a nearly zero denominator — the error this whole
paper is about.

### 2.3 Fieller as the two-by-one case

**Proposition.** Let `n = 2`, `k = 1`, `S` nonsingular, and let the model be the ray `v = (1, R₀)`. Then

```text
D(v) = z_F²,
```

where `z_F` is Fieller's statistic for the null `E[y₂]/E[y₁] = R₀`, namely

```text
z_F = (y₂ − R₀ y₁) / sqrt( S₂₂ − 2 R₀ S₁₂ + R₀² S₁₁ ).
```

*Proof.* Write `c = (−R₀, 1)ᵀ`, so `cᵀ v = 0` and `cᵀ y = y₂ − R₀ y₁`. Since `n − k = 1`, the orthogonal
complement of `span(v)` in the `S⁻¹` geometry is one-dimensional, and the minimized quadratic form equals the
squared `S`-standardized component of `y` along the annihilator of `v`:

```text
D(v) = (cᵀ y)² / (cᵀ S c).
```

The denominator is `cᵀ S c = R₀² S₁₁ − 2 R₀ S₁₂ + S₂₂`, which is Fieller's variance, and the numerator is
Fieller's contrast squared. ∎

The proposition is elementary, and we state it because it settles what the projective statistic *is*. It is not a
rival to Fieller's theorem; on Fieller's own problem it is Fieller's theorem, and it inherits Fieller's
guarantees, including the possibly-unbounded confidence set that Gleser and Hwang [3] show is unavoidable. What
it adds is that the same expression continues to make sense when `n > 2` and `k > 1`, where Fieller's contrast
has no analogue, because there is no longer a single `R₀` to invert over.

`tests/test_projective_inference.py` verifies the identity numerically at `R₀ ∈ {1, 4, 10.99, 16, 120.8,
2080.3}` — six orders of magnitude — to a relative agreement of `10⁻¹⁴`, and §4.3 verifies it on the real data.

### 2.4 Confidence regions without a chosen denominator

Inverting `D` over the model gives a confidence region directly in the space of directions. For a
one-parameter family of rays `v(θ)`, the set `{θ : D(v(θ)) ≤ χ²_{ν,1−α}}` is a confidence region for `θ` that
never requires a ratio to be formed and never requires a coordinate to be nonzero. Where the family is
`v(R₀) = (1, R₀)` this is Fieller's interval and can be unbounded; where the family is a genuinely
higher-dimensional model space, the region lives on a Grassmannian and unboundedness is replaced by the region
covering a positive-measure set of directions. The reporting discipline is the same in both cases: report the
region, not a point estimate with a symmetric error bar.

### 2.5 Annihilating functionals, and why they are the safest thing in the paper

There is a special case worth isolating because it removes not only the denominator but the matrix inverse.

Suppose a *class* of competing models — not one model, a class — all lie in some subspace `W ⊂ Rⁿ`. Then any
linear functional `c` with `c ⊥ W` annihilates the entire class, and

```text
cᵀ y ± sqrt(cᵀ S c)
```

is a complete test of the class, with no minimization, no pseudo-inverse, no chosen denominator, and no
degradation whatsoever when one rung sits close to zero. The variance of `cᵀ y` is a fixed quadratic form in `S`;
nothing about it becomes ill-conditioned.

For a ladder on rungs `r = 1, 2, 4`, the divided difference

```text
f[1,2,4] = ( m(4) − 3 m(2) + 2 m(1) ) / 6
```

is such a functional. It returns exactly `1` on `m(r) = r²`, exactly `0` on `m(r) = α + βr` for every `α, β`, and
by linearity it measures the coefficient of the quadratic part of any second-degree fit. It is the discrete
curvature of the ladder. Its sign is a class-level statement: a family of laws that is linear in `r` predicts
zero, a family that is convex predicts positive, and a *negative* measured value is incompatible with every
member of both classes at once, whatever their amplitudes.

We stress the epistemic difference between §2.2 and §2.5. The `χ²` statistic answers "is this particular model
excluded"; the annihilating functional answers "can the data have come from anything in this class at all". The
second question survives a much wider range of things going wrong with the first, and it is the one that carries
the physics in §4.5.

---

## 3. Why the ratio framework fails in a way that is hard to notice

### 3.1 Fieller fixes the distribution; it does not fix the conditioning

Replacing `(R̂ − R₀)/SE(R̂)` with Fieller's contrast repairs the coverage. It does not make the denominator
strong. In the application of §4, the nominated denominator sits `3.62σ` from zero (Table T1), which is enough
for the contrast to be well behaved but not enough for the inverted interval to be informative: the 3σ Fieller
interval for the deciding ratio is `[2.395, 27.468]`, of width 25.1 (Table T1). That interval contains four of
the eight competitors, including two whose rays the three-rung test excludes above 7σ.

The interval is not wrong. It is the honest interval for the ratio, given a denominator resolved to 3.6σ. The
error is upstream: the ratio should not have been the object of inference.

### 3.2 One denominator means one ratio at a time

A ladder of `n` rungs offers `n − 1` independent ratios against a nominated denominator. It is tempting to test
all of them. But they share the denominator, so they are strongly dependent, and their joint law is not described
by the marginal Fieller intervals. Combining them correctly requires the joint distribution of the contrast
*vector* `(y₂ − R₀,₂ y₁, …, y_n − R₀,ₙ y₁)` in the covariance of `y` — at which point one has written down
`D(v)` with `v = (1, R₀,₂, …, R₀,ₙ)` and abandoned the ratio framework anyway.

In practice the combination step is skipped, and a single ratio is nominated as *the* deciding entry. This is
the second discretionary choice, after the choice of denominator, and it is the one that governs §3.3.

### 3.3 The self-sealing move

Which ratio should decide? The analyst who knows the ladder knows which rungs are contaminated, and will
nominate the deciding entry to avoid them. In §4 the design did precisely this, for a good and explicitly stated
reason: the middle rung carries a geometric leakage from an unwanted angular harmonic, of opposite sign to the
leakage on the other rungs, and the design's own frozen text records the leakage coefficient exactly
(`1148/21025 = 0.0546`) and the argument for treating it as negligible.

The argument, traced to its source, is a model-selection result on a different observable — a homology character
— quoted as bounding the *angular* amplitude ratio `|A₈/A₄|` "well below 1". That inference is a plausibility
argument. It is not quantified anywhere, and it was never measured.

The self-sealing is now visible. The rung that was excluded from the decision is the only rung that carries
information about the harmonic whose smallness justified excluding it. A framework that must nominate one ratio
has no way to say "this rung disagrees with the assumption under which I dropped it," because the rung is not in
the statistic. A framework that tests a ray keeps all rungs in the statistic by construction, and the
disagreement appears as a residual.

---

## 4. Application: a three-rung aspect-ratio ladder at 580 sites

### 4.1 The measurement

The observable is the spin-4 amplitude `A₄` of a matching-odd percolation readout on a torus, extracted at three
aspect ratios `r = 1, 2, 4` at fixed site count 580, with 100 aligned batches per rung sharing a seed and replica
offset. Full design, engine and provenance are in `predictions/aspect_ladder_n580_20260905.yaml` and
`results/aspect-ladder-n580/latest.json`; the covariance is a paired delete-one jackknife in which batch `b` is
deleted simultaneously across all three rungs, which is valid precisely because the rungs share the stream.

The measured response and its per-rung resolution are Table **T1**:

```text
r = 1:  9.016433e-04 ± 2.4914e-04    3.62 σ from zero
r = 2:  2.910977e-03 ± 2.5911e-04   11.23 σ
r = 4:  4.131808e-03 ± 1.9513e-04   21.17 σ
```

The rung the design nominated as denominator is the weakest of the three by a factor of six.

### 4.2 The eight competitors

Each competitor is a ray `v = (m(1), m(2), m(4))` fixed before the reveal, in the frozen prediction file (T3):
two families linear in `r` (`bare_aspect_ratio` `(1,2,4)` and `no_modulus_dependence` `(1,1,1)`), plain area
scaling `(1,4,16)`, the weight-4 modular prediction `(1, 2.75, 10.9908)`, and weight-8 and weight-12 modular
laws reaching `(1, 32.5, 2080.3)`.

### 4.3 The two-entry control: the statistic did not change

Before comparing verdicts we must rule out the trivial explanation that the two procedures differ because they
are different statistics. Table **T2** restricts the projective statistic to exactly the two rungs the frozen
test used — `r = 1` and `r = 4`, the design's declared deciding entry — and compares it with Fieller's `z`
squared on the same pair, competitor by competitor, on the real covariance:

```text
largest relative deviation over all eight competitors:  1.5e-15
```

The two are the same number. Every verdict change reported below is therefore attributable to the third rung and
to nothing else.

### 4.4 The verdicts

Table **T3**. Frozen (one entry, Fieller) against projective (three rungs):

| competitor | frozen σ | projective σ | frozen | projective |
|---|---:|---:|---|---|
| `bare_aspect_ratio` | 0.50 | 2.74 | compatible | compatible |
| `q4_jordan_weight4` | 2.08 | **7.00** | compatible | **excluded** |
| `plain_area_scaling` | 2.56 | **7.13** | compatible | **excluded** |
| `weight8_E8` | 3.48 | 10.47 | excluded | excluded |
| `weight12_E4_cubed` | 3.61 | 11.25 | excluded | excluded |
| `weight12_E12` | 3.61 | 11.25 | excluded | excluded |
| `no_modulus_dependence` | 9.48 | 9.22 | excluded | excluded |
| `weight12_delta` | 21.17 | ∞ | excluded | excluded |

Two verdicts flip from compatible to excluded, at 7.0σ and 7.1σ. This is not a marginal reclassification: the
frozen scoring reported the weight-4 modular law as the best-supported surviving competitor at 2.08σ, and the
three-rung test excludes it above 7σ.

It matters for the credibility of the correction that it does not move everything in one direction. The
exclusion of `no_modulus_dependence` *weakens*, 9.48σ to 9.22σ, because the third rung is slightly more
consistent with a flat law than the deciding pair alone was. And `bare_aspect_ratio`, which the frozen test
reported as an essentially perfect fit at 0.50σ, becomes a 2.74σ tension. A correction that improved every number
in the direction of the author's preferred conclusion would deserve suspicion; this one does not have that
shape.

### 4.5 The curvature: nothing in the class can fit

Table **T4**. The annihilating functional of §2.5, applied to the same response:

```text
f[1,2,4]  =  −4.6631e−04  ±  1.5297e−04        z = −3.05
```

and, across every positive-definite value of the one covariance entry the source artifact did not store
(§5.2), `z ∈ [−3.79, −2.62]`. The sign never flips.

What the competitors predict for the same functional:

```text
bare_aspect_ratio      0 (exactly)      linear in r
no_modulus_dependence  0 (exactly)      linear in r
weight12_delta        +0.271
q4_jordan_weight4     +0.790
plain_area_scaling    +1.000
weight8_E8           +16.685
weight12_E4_cubed   +211.212
weight12_E12        +330.793
```

Every competitor predicts zero or strictly positive. The measurement is negative at 3.05σ, and by §2.5 this is a
statement about the *classes*, not about eight particular amplitudes: the response is concave in `r`, and neither
the linear class nor the convex modular class contains a concave member at any amplitude.

This is the paper's substantive physical finding, and it is the one statement here that is completely free of
denominators, matrix inverses and conditioning arguments. It is a single linear contrast with a fixed quadratic
form for its variance.

### 4.6 The design's assumption, priced

Section 3.3 described the move that dropped the middle rung: a spin-8 leakage of coefficient `0.0546`, assumed
small because `|A₈/A₄| ≪ 1`. We can now price that assumption. Fixing each competitor's amplitude from the two
clean rungs — whose leakages share a sign — and asking what spin-8 amplitude is needed to bring the middle rung
onto the predicted ray gives Table **T5**:

| competitor | required \|A₈/A₄\| |
|---|---:|
| `no_modulus_dependence` | **0.4** |
| `bare_aspect_ratio` | 7.7 |
| `plain_area_scaling` | 32.0 |
| `q4_jordan_weight4` | 32.1 |
| `weight8_E8` | 182 |
| `weight12_delta` | 222 |
| `weight12_E4_cubed` | 783 |
| `weight12_E12` | 785 |

Seven of the eight require the assumed-negligible harmonic to be between 8 and 785 times the spin-4 amplitude it
was assumed to be far smaller than. The only competitor consistent with a small ratio, `no_modulus_dependence` at
0.4, is excluded at 9.2σ by the two clean rungs alone — the one exclusion that does not depend on the disputed
rung at all.

So the ladder presents a fork with no third branch:

- **either** the spin-8 amplitude is one to three orders of magnitude larger than the bound under which the
  middle rung was set aside, in which case the frozen design's justification for setting it aside is false;
- **or** the assumed decomposition `C + A₄cos4θ + A₈cos8θ` is not the right form for this observable, in which
  case the amplitudes being compared are not the amplitudes the models predict.

The ratio framework could not present this fork, because the quantity that distinguishes the branches lives in
the residual direction the ratio removes.

---

## 5. What the reanalysis does and does not establish

### 5.1 Established

- On this data, the projective and Fieller statistics agree exactly where they are both defined (`1.5·10⁻¹⁵`).
- Using all three rungs, seven of eight prospectively frozen competitors are excluded at 7.0–11.3σ.
- The response's second divided difference is negative at 3.05σ, and no competitor's class permits that sign.
- Reconciling the middle rung with any competitor requires `|A₈/A₄|` between 0.4 and 785, and the only value
  below 1 belongs to a competitor excluded at 9.2σ without that rung.

### 5.2 One verdict is not stable

The source artifact was written by a scorer that stored the two pairwise covariances the ratios needed and not
the full `3×3` matrix; `cov(r₂, r₄)` is therefore absent. We do not assume it is zero — assuming zero is a
choice, and an undeclared one. Instead we scan every value that keeps the matrix positive definite,
`corr(r₂, r₄) ∈ [−0.985, 0.980]`, and report the range (T3):

```text
q4_jordan_weight4     σ ∈ [5.93, 8.86]     excluded throughout
plain_area_scaling    σ ∈ [6.10, 8.93]     excluded throughout
weight8_E8            σ ∈ [10.03, 10.97]   excluded throughout
weight12_*            σ ∈ [11.13, 11.38]   excluded throughout
no_modulus_dependence σ ∈ [9.21, ∞]        excluded throughout
bare_aspect_ratio     σ ∈ [2.14, 4.77]     NOT stable
```

Exactly one verdict — `bare_aspect_ratio`, the only surviving competitor — depends on the missing entry. That
verdict is therefore reported as undetermined, and it is worth exactly one deterministic replay under the
current scorer and nothing more. The scorer has since been changed to persist the full response vector and
covariance (§7), so this gap does not recur.

### 5.3 Not established

- **No competitor is supported.** Seven are excluded and the eighth is undetermined. "The data reject every law
  we froze" is the finding; it is not evidence for a ninth law, and in particular is not evidence for the
  surviving competitor, which fits three points and is a law that fits three points.
- **The two branches of §4.6 are not separated.** Two orientations per rung determine `C` and `A₄` exactly with
  nothing left over, so this measurement has no residual with which to check the angular form itself. A large
  spin-8 amplitude and a wrong functional form produce the same arithmetic here.
- **No modular identification.** Nothing here bears on whether the underlying object carries modular weight; the
  competitors are rays, and rays are compatible with many mechanisms.
- **No claim about the percolation threshold.** This is an amplitude ladder, not a threshold estimate.

---

## 6. The design consequence

A negative result that names its own next measurement is worth more than one that does not, and §4.6 names one
precisely.

The obstruction is that two lattice orientations per family determine two unknowns, `C` and `A₄`, with zero
degrees of freedom left over. The remedy is three orientations per family, which determine `C`, `A₄` and `A₈`
together and turn the assumed bound into a measured number. The constraint is arithmetic rather than
computational: an orientation of a torus family corresponds to a representation of a fixed norm as a Gaussian
integer, and three inequivalent representations are needed at *both* the square and the rectangular norm
simultaneously.

The smallest site count meeting that condition is **N = 650**: square `|w|² = 650` with `25+5i`, `23+11i`,
`19+17i`, and rectangular `|w|² = 325` with `18+i`, `17+6i`, `15+10i`. This is 2.24× the sites of the design in
which the leakage was first identified, and requires three coupled runs rather than two. It is the measurement
this reanalysis makes mandatory rather than optional, and it converts the fork of §4.6 into a decision: a
measured `|A₈/A₄|` near 0.05 falsifies the large-leakage branch and indicts the angular form; a measured value
above 1 falsifies the design's own justification for the analysis it froze.

---

## 7. Reproducibility and the estimator's preconditions

### 7.1 Software

| component | file |
|---|---|
| projective statistic, pseudo-inverse, `χ²` tail | `scripts/projective_inference.py` |
| the N=580 rescoring | `scripts/aspect_ladder_projective_rescore.py` |
| this manuscript's tables | `scripts/p3_manuscript_evidence_table.py` |
| Fieller↔projective identity, 8 tests | `tests/test_projective_inference.py` |
| rescoring, 11 tests | `tests/test_aspect_ladder_projective_rescore.py` |
| table generation | `tests/test_p3_manuscript_evidence_table.py` |

The implementation uses `mpmath` at 50 digits and no other numerical dependency. Every table in this draft is
regenerated from `results/p3-projective-inference-manuscript/latest.json` and a regression test fails if the
rendered tables drift from the artifact.

### 7.2 What the covariance estimate requires

The method's one real precondition is that `S` estimate the covariance of a single random object across all
coordinates. Three practices make that true and are cheap:

1. **Share the random stream across rungs.** All rungs use one seed and one replica offset, so a delete-one
   jackknife can delete batch `b` from every rung simultaneously and the resulting pseudo-values are paired. A
   ladder whose rungs were run independently has a diagonal `S` by construction, which is not a measurement of
   independence but an absence of one.
2. **Persist the response vector and its full covariance, not the derived comparisons.** The gap of §5.2 exists
   only because the original scorer stored what its ratios needed. Storing `y` and `S` costs `O(n²)` numbers and
   makes every future reanalysis possible without a rerun; storing derived ratios does not. This is the concrete
   form, for amplitude ladders, of a general rule: **do not quotient the amplitude away in storage.**
3. **Report the retained rank of `S` with every `D`.** See §2.2.

### 7.3 Chronology

The design, its competitor list, its sample counts, its stopping rule and its declared systematic were frozen in
`predictions/aspect_ladder_n580_20260905.yaml` before the run; the reveal, the Fieller correction and the
covariance replay are separately dated in the repository's evidence ledger; the reanalysis of this paper adds no
samples and rescores a committed artifact.

---

## 8. Relation to existing work

Fieller's theorem [1] is the correct treatment of a single ratio and, as §2.3 shows, the two-dimensional case of
the statistic advocated here. Hinkley [2] and Marsaglia [4] give the exact distribution of a ratio of correlated
normals and make quantitative the tail behaviour that motivates avoiding it. Gleser and Hwang [3] establish that
the unboundedness of Fieller's region is unavoidable rather than an artifact, which is the theoretical reason a
finite symmetric error bar on a ratio should be treated as a red flag. Quadratic forms in normal variables with
generalized inverses, including the rank-based degrees of freedom used in §2.2, are standard [5].

What we have not found in the finite-size scaling literature is the routine use of the projective form. The
practice of testing amplitude ratios against predicted ratios is ubiquitous — it is how conformal-invariance
predictions for aspect-ratio and boundary-condition dependence are compared with lattice data [6] **[LIT]** —
and where the denominator is strong nothing goes wrong. The contribution here is to identify the regime where
something does go wrong, to show that the failure is a design failure rather than an arithmetic one, and to give
a drop-in replacement that specializes to the familiar procedure whenever the familiar procedure was valid.

The annihilating-functional idea of §2.5 is a special case of testing a linear hypothesis on a linear model, but
we have not seen it deployed as a *class-level* tool for finite-size ladders, where its virtue is that it
requires no model to be nominated and no matrix to be inverted. Divided differences are the natural annihilators
for polynomial classes in the rung variable; other classes admit others, and constructing them for a given
competitor family is a small exercise in linear algebra rather than a new method.

---

## References

1. Fieller, E. C. (1954). "Some problems in interval estimation." *Journal of the Royal Statistical Society B*
   **16**, 175–185.
2. Hinkley, D. V. (1969). "On the ratio of two correlated normal random variables." *Biometrika* **56**, 635–639.
3. Gleser, L. J. and Hwang, J. T. (1987). "The nonexistence of 100(1−α)% confidence sets of finite expected
   diameter in errors-in-variables and related models." *Annals of Statistics* **15**, 1351–1362.
4. Marsaglia, G. (2006). "Ratios of normal variables." *Journal of Statistical Software* **16**, 1–10.
5. Rao, C. R. and Mitra, S. K. (1971). *Generalized Inverse of Matrices and its Applications*. Wiley.
6. **[LIT]** Cardy, J. L. (1984–1986), finite-size scaling and conformal invariance; specific citation to be
   fixed against the primary when the aspect-ratio amplitude comparison of §8 is stated in the final draft.
7. Quenouille, M. H. (1956). "Notes on bias in estimation." *Biometrika* **43**, 353–360; Tukey, J. W. (1958).
   "Bias and confidence in not quite large samples" (abstract). *Annals of Mathematical Statistics* **29**, 614.
8. Efron, B. and Stein, C. (1981). "The jackknife estimate of variance." *Annals of Statistics* **9**, 586–596.
