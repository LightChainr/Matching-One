# Research Synthesis — 2026-08-28

This is the execution-facing synthesis of Matching One. It states what is strongest, what just failed, and which experiments now have the highest information value.

## Current thesis

The project is no longer best described as one pure two-spin-4 finite-size model.

The evidence now separates into three layers:

```text
matching-odd central orientation sector
    -> robust and prospectively reproduced
    -> compatible with DeltaCos4 * N^-13/8
    -> local root movement through DeltaRoot ~= -DeltaM/M'

matching-even central sector
    -> original positive N^-1 assignment prospectively fails in sign

S-prime derivative sector
    -> definitely nonzero
    -> pure N^-5/4 fails
    -> finite-size correction required
```

The first line is currently the strongest scientific result in the repository. The second and third are now model-selection problems rather than support for the original conjunction.

## Strongest evidence

### 1. Independent same-N orientation confirmation

At `N=65,85,130,145,170`, a fresh 100M-per-size seed reproduces the orientation sign predicted by `Delta cos(4 theta)` at all five sizes:

```text
z = 16.03, 11.23, 5.22, 5.27, 2.58
A4 = N^(13/8) DeltaM/DeltaCos4 = 0.7885 +/- 0.0352
chi2 = 1.53 / 4
```

This establishes a reproducible finite-size orientation sector independently of any continuum interpretation.

### 2. Three prospective norm-2 Gaussian lineages

For raw orientation contrasts, the H4-like law predicts

```text
DeltaM(2N)/DeltaM(N) = -2^(-13/8) = -0.3242098887...
```

and the three exact Gaussian genealogies are compatible with that transformation:

```text
65 -> 130 = -0.31382 +/- 0.0908
85 -> 170 = -0.34095 +/- 0.1118
145 -> 290 fixed-p child residual z = -0.483
```

`1+i` cannot distinguish H4 from H12/H20, but it is a strong parameter-free test of the odd square-harmonic semigroup class.

### 3. Prospective new geometries N=185/265 preserve the odd law

The 500M-per-size full-curve run is the most important new result because these geometries were unused target data.

```text
N=185 DeltaM = +1.36830e-4 +/- 3.42352e-5
N=265 DeltaM = +1.27110e-4 +/- 3.46783e-5

x=21/4 H4: chi2 = 3.04598 / 2
zero:       chi2 = 29.40938 / 2
x=17/4:     chi2 = 30.24613 / 2
```

The frozen x=21/4 H4-like radial law therefore survives a prospective new-geometry challenge and clearly outperforms both zero and the larger x=17/4 adversary.

This materially strengthens the surviving central odd-sector claim. It still does not make H4 unique.

### 4. Root movement remains locally tied to DeltaM

Threshold-rank reconstructions give

```text
-DeltaRoot * mean(M') / DeltaM ~= 1
```

on the tested systems. P45 also passes a frozen angular-normalized root-amplitude test at N=65/85. The empirical residual therefore moves the finite matching root through the expected local mechanism.

### 5. Full-curve data resolve finite-size slope structure

Clean 100M norm-2 full curves show that the simple bare slope multiplier is too precise a statement:

```text
observed ratios = 1.2939835, 1.2943776
2^(3/8)         = 1.2968396
chi2            = 6412.89 / 2
```

The discrepancy is only about 0.2%, but it is real at this precision. The natural interpretation is a finite-size thermal-metric correction, not abandonment of the leading `y_t=3/4` scaling law.

A pre-target minimal relative-`1/N` correction is already frozen for the third `145->290` full-curve lineage under #50.

## The important failure: the simple even sector reverses sign

The same prospective N=185/265 run falsifies the original positive matching-even central prediction:

```text
N=185 DeltaS = -6.08154e-5 +/- 8.08957e-6
N=265 DeltaS = -7.02495e-5 +/- 9.38562e-6

frozen positive N^-1: chi2 = 240.24721 / 2
zero:                  chi2 = 112.53891 / 2
```

Both observed effects are about `-7.5` sampling SE from zero. This is not a small amplitude miss and should not be explained away as noise.

Consequences:

- do not describe the original two-sector model as confirmed;
- do not attach the even-sector failure to the surviving odd `DeltaM` law;
- first re-audit the signed `P4[S]` sequence with the corrected normalized-P4 convention;
- use existing data before commissioning another run.

Issue #48 now owns this question.

## The derivative result: S-prime is real, but pure scaling fails

On the same new geometries,

```text
Y = N^(5/4) P4[S']
N=185: 2.57971
N=265: 2.85844
```

and the frozen scores are

```text
pure N^-5/4:  chi2 = 52.71634 / 2
rank-2/log:    chi2 =  1.20360 / 2
analytic 1/N:  chi2 =  0.86221 / 2
zero:          chi2 = 1278.55524 / 2
```

So `S'` is decisively nonzero, the pure law is prospectively rejected, and both fixed correction forms survive. The data do not yet select a unique q=2 versus logarithmic/Jordan mechanism.

This is an example of why the project should retain chronological fixed-model scores instead of choosing whichever retrospective curve looks best.

## What remains genuinely weak

### H4 versus H12

`1+i` multiplication only identifies the odd square-harmonic class. H4, H12, H20, ... all reverse under a pi/4 rotation.

Norm 5 breaks that alias. For the frozen N=325/425 children, raw H4 and H12 predict opposite signs. This is now the single highest-information operator-level experiment.

### Unique asymptotic radial exponent

The `13/8` law has become substantially harder to dismiss: it now survives independent seeds, old holdouts, three exact norm-2 lineages, and new N=185/265 geometries while beating the x=17/4 adversary there.

But the synthetic red-team still warns that finite sets of N have weak power against some corrected radial mechanisms. New multiplier leverage remains better than just shrinking old error bars.

### Continuum operator identification

The `x=21/4`, spin-4 thermal-family candidate has the right arithmetic and remains the leading interpretation of the odd central sector, but unique harmonic content, matching/OPE parity, logarithmic structure, and lattice coupling are not proven.

Empirical matching-even/odd decomposition is weaker than a full local OPE/interchiral automorphism. Keep those claims separate.

## If we can run only two expensive next experiments

### 1. Norm-5 Gaussian spectroscopy — #57

Highest priority.

The production engine already supports the frozen children

```text
N=325: (17,6) - (18,1)
N=425: (16,13) - (19,8)
```

Raw child/parent predictions include

```text
H4  = -0.04096017184...
H8  = -0.12334863177...
H12 = +0.11003540563...
```

H4 and H12 predict opposite child signs.

Start with a 1M threshold-rank variance pilot, extend the pilot only if the SE estimate is unstable, and choose final production size from measured discrimination power. Do not assume multi-billion replicas in advance.

Reuse the same full curves for normalized derivatives, root gaps, and q=2-versus-log transfer.

### 2. Third full-curve lineage — #50

Score `145 -> 290` on the full curve.

Order:

1. raw residual law;
2. bare slope baseline;
3. already-frozen finite-size slope correction;
4. raw and corrected induced root targets;
5. derivative channels with the normalized-P4 sign convention.

This tests whether the finite-size correction resolved independently in the first two norm-2 lineages predicts a third lineage.

## Zero-extra-compute work now worth doing

### Even/derivative sequence — #48

Before any dedicated new production:

- recompute signed `P4[S], P4[D], P4[D'], P4[S']` across N=65..265 directly from the stored full curves;
- verify that the prospective S sign reversal is not an orientation-order bookkeeping artifact;
- examine simple crossing/correction forms for S;
- retain q=2 and log/Jordan as the first fixed S-prime correction models;
- use cross-channel covariance to test whether S and S-prime corrections share a finite-size source.

Then reuse #50/#57 as fresh leverage.

## Secondary tracks

The following are useful, but should not displace #57/#50:

- axis-annihilator q=3/V14 test if its exact-marginal CRN engine/scorer remains clean;
- bounded exact matching-polynomial complex-zero work;
- exact N=1105 H0/H4/H8/H12 projector, gated behind cheaper multipliers;
- kappa3/continuum bridge;
- literature completion;
- bounded PSLQ only after provenance constraints.

## What to stop doing

Do not spend major effort on:

- more N=185/265 replicas merely to rescue the failed even conjunction;
- another five-size free-exponent fit;
- theory notes without a distinct frozen prediction or exact test;
- broad PSLQ searches;
- N=1105 before norm 5;
- generalized infrastructure rewrites that do not unlock a current experiment.

## Current project thesis

> Square-site/matching finite-size corrections contain a robust matching-odd orientation sector. Its central residual is compatible with `DeltaCos4*N^-13/8` across independent seeds, exact Gaussian transformations, and prospective new geometries, and it moves the finite root through the expected local mechanism. The original simple matching-even companion law is prospectively falsified in sign, while the S-prime derivative channel is nonzero but requires finite-size correction. The immediate mathematical task is therefore to identify the odd harmonic with norm-5 spectroscopy and to resolve the even/derivative correction structure without weakening the successful odd-sector evidence.

That is the line to optimize around until #57 or #50 breaks it.
