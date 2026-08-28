# Research Synthesis — 2026-08-28

This is the execution-facing synthesis of Matching One after the N=185/265 prospective reveal and the subsequent cross/either channel-map erratum.

## Current thesis

The strongest current finite-size picture is a **two-parity orientation structure with one specifically corrected derivative channel**.

```text
matching-odd central sector
    -> robust and prospectively reproduced
    -> compatible with DeltaCos4 * N^-13/8
    -> local root movement through DeltaRoot ~= -DeltaM/M'

matching-even central sector
    -> compatible with the frozen N^-1 amplitude after exact either->cross channel conversion

intrinsic-center P48 spectrum on new geometries
    P4[S]  ~ N^-1      -> survives
    P4[D]  ~ N^-13/8   -> survives
    P4[D'] ~ N^-5/8    -> survives
    P4[S'] ~ N^-5/4    -> pure law fails; correction required
```

The main scientific uncertainty is now operator/correction identification, not existence of the leading orientation signal or broad empirical parity structure.

## Strongest evidence

### 1. Same-N orientation signal and Gaussian semigroup

Independent high-statistics N=65,85,130,145,170 data reproduce the signed orientation signal and are compatible with

```text
DeltaM ~ DeltaCos4 * N^-13/8.
```

Three prospective `1+i` Gaussian genealogies are compatible with the raw no-fit transformation

```text
DeltaM(2N)/DeltaM(N) = -2^(-13/8).
```

The third fixed-coordinate 145->290 child was frozen before its run and passed with residual z about -0.48.

### 2. N=185/265 prospective new geometries support the odd radial law

The 500M-per-size target block gives

```text
N=185 DeltaM = +1.36830e-4 +/- 3.42352e-5
N=265 DeltaM = +1.27110e-4 +/- 3.46783e-5

x=21/4 H4-like: chi2 = 3.04598 / 2
zero:            chi2 = 29.40938 / 2
x=17/4:          chi2 = 30.24613 / 2
```

This is the strongest genuinely new-geometry support for the current matching-odd law. H4 itself is still not unique because odd higher harmonics remain possible until norm-5 spectroscopy.

### 3. Matching-even central N=185/265 agrees after the exact channel map

The initially published #108 interpretation said the even sector reversed sign. A protocol audit found that the frozen source and target statistics were different wrapping channels:

- frozen amplitude: P31 `either/even`;
- target threshold rank: `cross/even`.

P31 already showed these channels have opposite orientation contrasts, and complementary torus topology gives

```text
DeltaS_cross = -DeltaS_either.
```

The channel-corrected frozen predictions are

```text
N=185: -6.75216374588e-5
N=265: -6.89194469703e-5
```

against observations

```text
N=185: -6.08153762334e-5 +/- 8.08956556558e-6
N=265: -7.02495078452e-5 +/- 9.38562007761e-6.
```

With the unchanged source covariance,

```text
corrected chi2 = 0.5700315436 / 2
z = +0.667, -0.119.
```

No target fit is introduced. PR #134 adds the exact regression and erratum while preserving the original #108 artifacts.

This restores empirical consistency of the matching-even `N^-1` amplitude over the new geometries. It does **not** prove the identity-family `x=4` operator assignment or a local OPE automorphism.

### 4. Prospective intrinsic-center P48 parity score

The same N=185/265 full curves are independent target geometries for the P48 amplitudes frozen from N=65,85,130. Scoring with zero target refits gives

```text
P4[S]   ~ N^-1:     chi2 =  1.13878 / 2   zero = 112.53974 / 2
P4[D]   ~ N^-13/8:  chi2 =  0.28085 / 2   zero =  29.40844 / 2
P4[D']  ~ N^-5/8:   chi2 =  0.08761 / 2   zero =  59.39319 / 2
P4[S']  ~ N^-5/4:   chi2 = 52.71634 / 2   zero = 1278.55524 / 2
```

Thus `S`, `D`, and `D'` pure laws all transfer to genuinely new geometries. `S'` is the unique clear pure-law failure among these four intrinsic-center channels.

This is distinct from the Issue #43 fixed-coordinate cross/either erratum: the P48 score uses intrinsic-center projectors and their own frozen amplitudes. Together, the two analyses strengthen the empirical parity picture rather than merely repairing a sign mistake.

### 5. Root movement remains tied to the measured residual

Threshold-rank reconstructions give

```text
-DeltaRoot * mean(M') / DeltaM ~= 1
```

on the tested systems, and P45 passes an angular-normalized root-amplitude test. The finite matching root therefore moves through the expected local residual/slope mechanism.

### 6. The derivative sector exposes real subleading physics/statistics

Clean 100M full curves resolve that the bare center-slope multiplier is not exactly `2^(3/8)` at current N. The discrepancy is only about 0.2% but statistically decisive, so finite-size thermal-metric corrections are measurable.

For the one prospective pure-law failure, `S'`, the N=185/265 scores are

```text
pure P4[S'] ~ N^-5/4: chi2 = 52.71634 / 2
rank-2/Jordan log:      chi2 =  1.20360 / 2
analytic 1/N:           chi2 =  0.86221 / 2
zero:                    chi2 = 1278.55524 / 2
```

Thus `S'` is real, its pure leading power is insufficient, and both predeclared correction mechanisms survive. This remains a genuine unresolved problem after the channel erratum.

## What is still weak

### H4 versus H12/H20

Norm-2 `1+i` rotation identifies the odd square-harmonic class, not H4 uniquely. Norm-5 N=325/425 predicts different magnitudes and, for raw H4 versus H12, opposite signs. This is the most important expensive next test.

### Unique radial/operator identity

The `13/8` law is now difficult to dismiss over the tested range, but finite-size corrected radial alternatives must be tested with multiplier/new-size leverage, not another free fit on the same five points.

### q=2 versus Jordan correction

The new N=185/265 data rule out the pure `S'` law but do not distinguish the two frozen corrected forms. Further information should come from thermal-coordinate shape, new multipliers/moduli, or joint-channel constraints—not model proliferation.

### Continuum meaning of matching parity

Empirical pair-exchange sectors are established more strongly than any local RG/OPE automorphism. A direct FK/Potts, pivotal/four-arm, or exact topological derivation remains missing.

## Highest-value next work

### 1. Norm-5 Gaussian spectroscopy — #57

The engine supports the frozen children

```text
N=325: (17,6) - (18,1)
N=425: (16,13) - (19,8).
```

Use raw and normalized conventions explicitly. Raw H4 and H12 predictions differ in sign; normalized P4 removes the H4 angular factor and therefore uses positive `Q^-alpha` for pure H4. Do not mix these conventions.

Use a frozen variance/power pilot and consider a predeclared optional-stopping-safe likelihood/e-value rule for any billion-replica production.

### 2. Third full-curve lineage — #50

Score 145->290 in the frozen order:

1. central residual transfer;
2. bare slope baseline;
3. already-frozen finite-size slope correction;
4. induced root prediction;
5. derivative channels.

This determines whether the small but decisive P49 slope drift has a predictive correction structure.

### 3. Extract more structure from existing full curves

Highest-value zero/low-compute analyses:

- S-prime correction and fixed-coordinate/intrinsic-center mapping #48;
- prequential evidence ledger #95;
- pivotal/Russo bridge #100;
- intrinsic quantile-center spectroscopy #101;
- multi-u thermal-response discrimination #119;
- joint operator-mixing matrix #125;
- information-optimal Gaussian design #102.

These should reduce the live model space before asking for another large target run.

## Routes we were underusing

The post-P43 frontier program adds several orthogonal ways to test the mechanism:

- derive the observable through FK/Potts torus sectors (#114);
- connect the `3/4` slope and `21/4` candidate to pivotal/four-arm geometry (#100/#121);
- vary torus modulus, not only microscopic orientation (#103);
- build exactly-critical tunable anisotropy controls with isoradial/star-triangle models (#106);
- derive the matching identity at configuration level using Euler/Betti data (#111);
- search for universal amplitude ratios across microscopic controls (#118);
- falsify/extend the exact N=10 `Beta(3,3)` threshold law on N=26 (#115);
- certify finite-polynomial Galois complexity (#104).

Higher-risk routes—discrete holomorphic defects, correlated hyperedge self-duality, transfer-matrix eigenoperator spectroscopy, rigorous finite-size bounds, local complex zeros and full-distribution collapse—remain gated on cheap theory/control evidence.

## What to stop doing

Do not spend major effort on:

- interpreting a wrapping-channel sign without naming the channel;
- treating `S'` correction as failure of the entire P48 parity spectrum;
- adding another correction exponent after every residual;
- broad PSLQ searches;
- N=1105 before cheaper norm-5/new-geometry tests;
- transfer-matrix width solely for more decimal digits;
- GPU runs without a frozen model-discrimination target;
- treating several derived observables from the same histogram as independent evidence.

## Current project thesis

> Square-site/matching finite-size corrections contain reproducible matching-odd and matching-even orientation sectors once exact wrapping-channel semantics are aligned. The central matching-odd residual is compatible with an odd-square-harmonic `N^-13/8` law across independent seeds, exact Gaussian transformations and prospective new geometries, and it moves the finite root through the expected residual/slope mechanism. The matching-even central amplitude also survives the N=185/265 prospective block after exact either-to-cross conversion. Independently, the intrinsic-center P48 pure laws for `S`, `D`, and `D'` all transfer to the same new geometries; `S'` is the one clear pure-law failure and requires finite-size correction. The remaining hard questions are whether the leading harmonic is genuinely H4, whether the `x=21/4` thermal-family interpretation is correct, and which subleading mechanism explains the `S'` and slope corrections.

That is the line to optimize around until norm-5 or the third full-curve lineage breaks it.
