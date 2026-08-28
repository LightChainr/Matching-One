# Research Synthesis — 2026-08-28

This is the execution-facing synthesis of Matching One. It states what is strongest, what actually failed, and which experiments have the highest information value.

## Current thesis

The project now has a cleaner empirical decomposition than the first reading of the N=185/265 result suggested:

```text
central matching-odd orientation sector
    -> robust and prospectively reproduced
    -> compatible with DeltaCos4 * N^-13/8
    -> local root movement through DeltaRoot ~= -DeltaM/M'

intrinsic-center P48 parity channels
    P4[S]  ~ N^-1      -> prospectively survives
    P4[D]  ~ N^-13/8   -> prospectively survives
    P4[D'] ~ N^-5/8    -> prospectively survives
    P4[S'] ~ N^-5/4    -> pure law fails; correction required

Issue #43 registered fixed-coordinate DeltaS
    -> original artifact fails numerically
    -> but source=either/even and target=cross/even were mismatched
    -> exact channel-map repair gives chi2 0.570/2
    -> preserve as protocol correction, not retroactive preregistered pass
```

The immediate scientific problem is therefore **not** a broad failure of matching-even parity. It is narrower: identify the odd harmonic, understand the finite-size slope correction, and distinguish the correction mechanism in `S'` while keeping channel/thermal-coordinate contracts explicit.

## Strongest evidence

### 1. Independent same-N orientation confirmation

At `N=65,85,130,145,170`, a fresh 100M-per-size seed reproduces the orientation sign predicted by `Delta cos(4 theta)` at all five sizes:

```text
z = 16.03, 11.23, 5.22, 5.27, 2.58
A4 = N^(13/8) DeltaM/DeltaCos4 = 0.7885 +/- 0.0352
chi2 = 1.53 / 4
```

This establishes a reproducible finite-size orientation sector independently of continuum interpretation.

### 2. Three prospective norm-2 Gaussian lineages

For raw orientation contrasts,

```text
DeltaM(2N)/DeltaM(N) = -2^(-13/8) = -0.3242098887...
```

and the exact Gaussian genealogies are compatible with that transformation:

```text
65 -> 130 = -0.31382 +/- 0.0908
85 -> 170 = -0.34095 +/- 0.1118
145 -> 290 fixed-p child residual z = -0.483
```

`1+i` identifies the odd square-harmonic semigroup class but cannot distinguish H4 from H12/H20.

### 3. Prospective N=185/265 geometries preserve the odd radial law

The 500M-per-size full-curve run is the strongest new-geometry test:

```text
N=185 DeltaM = +1.36830e-4 +/- 3.42352e-5
N=265 DeltaM = +1.27110e-4 +/- 3.46783e-5

x=21/4 H4: chi2 = 3.04598 / 2
zero:       chi2 = 29.40938 / 2
x=17/4:     chi2 = 30.24613 / 2
```

The frozen x=21/4 H4-like law survives and clearly outperforms both zero and the tested x=17/4 adversary. This materially strengthens the central odd-sector claim; it still does not make H4 unique.

### 4. Prospective intrinsic-center P48 score: three pure laws survive

Using source amplitudes frozen from `N=65,85,130` and independent N=185/265 targets, with zero target refits:

```text
P4[S]   ~ N^-1:     chi2 =  1.13878 / 2
P4[D]   ~ N^-13/8:  chi2 =  0.28085 / 2
P4[D']  ~ N^-5/8:   chi2 =  0.08761 / 2
P4[S']  ~ N^-5/4:   chi2 = 52.71634 / 2
```

The first three are compatible with the frozen pure amplitudes and strongly beat zero. `S'` is the unique clear pure-law failure.

This is stronger evidence for the empirical parity decomposition than the old retrospective N=145/170 score alone.

### 5. S-prime is nonzero and requires correction

On the new geometries:

```text
Y = N^(5/4) P4[S']
N=185: 2.57971
N=265: 2.85844
```

Frozen correction scores:

```text
pure N^-5/4:  chi2 = 52.71634 / 2
rank-2/log:    chi2 =  1.20360 / 2
analytic 1/N:  chi2 =  0.86221 / 2
zero:          chi2 = 1278.55524 / 2
```

Both correction forms survive. The smaller descriptive analytic statistic does not identify the mechanism or reverse preregistered chronology.

### 6. Root movement remains locally tied to DeltaM

Threshold-rank reconstructions give

```text
-DeltaRoot * mean(M') / DeltaM ~= 1
```

on tested systems. P45 also passes a frozen angular-normalized root-amplitude test at N=65/85.

### 7. Full-curve data resolve finite-size slope structure

Clean 100M norm-2 full curves show the bare center-slope multiplier is not an exact finite-N law:

```text
observed ratios = 1.2939835, 1.2943776
2^(3/8)         = 1.2968396
chi2            = 6412.89 / 2
```

The roughly 0.2% discrepancy is real at this precision. A minimal relative-`1/N` correction is frozen prospectively for the third `145->290` full-curve lineage.

## Protocol correction: what the Issue #43 DeltaS failure really means

The pre-target Issue #43 artifact sourced its even amplitude from P31 **`either/even`**. The target threshold-rank engine and scorer reconstruct **`cross/even`**. The registered positive means therefore give the historical failure

```text
chi2 = 240.24721 / 2.
```

That failed preregistration artifact must remain unchanged.

However, P31 already contains both source channels and verifies the exact orientation-difference map

```text
DeltaS_cross = -DeltaS_either.
```

Applying only this pre-existing map, with zero target refits, gives

```text
corrected cross/even chi2 = 0.57003 / 2
marginal residual z = +0.667, -0.119.
```

This is a **post-reveal protocol correction**, not a retroactive preregistered pass. It changes the scientific interpretation: the `240/2` failure identifies a source/target channel-contract bug, not physical falsification of the matching-even cross N^-1 law.

The intrinsic-center P48 `P4[S]` score is a separate construction and prospectively passes.

## What remains genuinely weak

### H4 versus H12

`1+i` multiplication only identifies the odd square-harmonic class. Norm 5 breaks that alias: frozen N=325/425 raw H4 and H12 predictions have opposite signs. This is the highest-information operator-level experiment.

### Unique asymptotic radial exponent

The `13/8` law now survives independent seeds, old holdouts, three norm-2 lineages, and new N=185/265 geometries while beating the tested x=17/4 adversary. But finite-N corrected radial alternatives still deserve prospective leverage; new multipliers are more useful than simply shrinking old error bars.

### S-prime correction mechanism

The data require a correction, but do not yet choose ordinary q=2 versus logarithmic/Jordan structure. This is the cleanest unresolved derivative question.

### Continuum operator identification

The `x=21/4`, spin-4 thermal-family candidate remains the leading interpretation of the central odd sector, but unique harmonic content, matching/OPE parity, logarithmic structure, and microscopic lattice coupling are not proven.

## If we can run only two expensive next experiments

### 1. Norm-5 Gaussian spectroscopy — #57

Highest priority.

```text
N=325: (17,6) - (18,1)
N=425: (16,13) - (19,8)

raw H4  ratio = -0.04096017184...
raw H8  ratio = -0.12334863177...
raw H12 ratio = +0.11003540563...
```

H4 and H12 predict opposite child signs.

Start with a 1M threshold-rank variance pilot, extend only if the SE estimate is unstable, then choose final production size from measured discrimination power. Reuse the same full curves for normalized derivatives, root gaps, and S-prime q=2-versus-log transfer.

### 2. Third full-curve lineage — #50

Score `145 -> 290` in order:

1. raw residual semigroup law;
2. bare slope baseline;
3. frozen finite-size slope correction;
4. raw and corrected induced root targets;
5. P48 derivative channels under the canonical normalized-P4 convention.

This tests whether the finite-size correction learned from the first two clean norm-2 lineages predicts a third lineage.

## Zero-extra-compute work — #48

The target is now specific:

- keep the prospective four-channel P48 score machine-readable;
- preserve the original Issue #43 DeltaS failure and exact channel-map repair side by side;
- assemble the signed scaled `S'` sequence N=65..265 with covariance where available;
- retain q=2 and log/Jordan as first correction models;
- compare fixed `p_ref`, intrinsic-center, and derivative transport to understand channel/thermal-coordinate bookkeeping;
- reuse #50/#57 before commissioning any new geometry.

## Secondary tracks

Useful, but subordinate to #57/#50:

- axis-annihilator q=3/V14 program;
- bounded exact matching-polynomial complex-zero work;
- exact N=1105 H0/H4/H8/H12 projector;
- kappa3/continuum bridge;
- literature completion;
- bounded PSLQ.

## What to stop doing

Do not spend major effort on:

- more N=185/265 replicas;
- interpreting the old `240/2` DeltaS score as physical even-sector falsification;
- another five-size free-exponent fit;
- theory notes without a distinct frozen prediction or exact test;
- broad PSLQ searches;
- N=1105 before norm 5;
- generalized infrastructure rewrites that do not unlock an active experiment.

## Current project thesis

> Square-site/matching finite-size corrections contain a robust orientation sector whose central matching-odd residual is compatible with `DeltaCos4*N^-13/8` across independent seeds, exact Gaussian transformations, and prospective new geometries, and moves the finite root through the expected local mechanism. The intrinsic-center P48 parity pattern also survives prospectively in `S`, `D`, and `D'`; the specific unresolved channel is `S'`, which is nonzero but requires finite-size correction. The apparent Issue #43 even-sector sign failure was a registered `either/even` source versus `cross/even` target contract mismatch, not a physical falsification. The immediate tasks are norm-5 harmonic spectroscopy, a third full-curve lineage, and correction spectroscopy in `S'`.

That is the line to optimize around until #57 or #50 breaks it.
