# Research Synthesis — 2026-08-29 sync

This is the execution-facing synthesis of Matching One after the N=185/265 prospective reveal, the cross/either channel-map erratum, and the repository-organization/exact-control integration pass completed on 2026-08-29.

## Current thesis

The strongest finite-size picture remains a **two-parity orientation structure with a real derivative-correction problem**:

```text
matching-odd central sector
    -> robust and prospectively reproduced
    -> compatible with DeltaCos4 * N^-13/8
    -> local root movement through DeltaRoot ~= -DeltaM/M'

matching-even central sector
    -> compatible with the frozen N^-1 amplitude after exact either->cross channel conversion

intrinsic derivative sectors
    -> S, D and D-prime pure laws survive N=185/265
    -> S-prime is decisively nonzero
    -> pure N^-5/4 fails prospectively
    -> q=2 and Jordan-log corrections both remain live
```

The main scientific uncertainty is operator/correction identification, not existence of the leading orientation signal.

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

```text
N=185 DeltaM = +1.36830e-4 +/- 3.42352e-5
N=265 DeltaM = +1.27110e-4 +/- 3.46783e-5

x=21/4 H4-like: chi2 = 3.04598 / 2
zero:            chi2 = 29.40938 / 2
x=17/4:          chi2 = 30.24613 / 2
```

This is the strongest genuinely new-geometry support for the current matching-odd law. H4 itself is still not unique because odd higher harmonics remain possible until norm-5 spectroscopy.

### 3. Matching-even N=185/265 agrees after the exact channel map

The original #108 interpretation mixed a frozen P31 `either/even` source with a threshold-rank `cross/even` target. Complementary torus topology gives

```text
DeltaS_cross = -DeltaS_either.
```

With no target refit:

```text
corrected chi2 = 0.5700315436 / 2
z = +0.667, -0.119.
```

The original unconverted score remains historical provenance. The correction is now enforced by the typed observable/channel layer, so a scorer cannot silently compare mismatched channels.

### 4. The prospective intrinsic-center four-channel picture is now reproducible from a canonical scorer

Using only frozen N=65/85/130 source amplitudes and independent N=185/265 targets:

```text
P4[S]   ~ N^-1:     chi2 =  1.13878 / 2
P4[D]   ~ N^-13/8:  chi2 =  0.28085 / 2
P4[D']  ~ N^-5/8:   chi2 =  0.08761 / 2
P4[S']  ~ N^-5/4:  chi2 = 52.71634 / 2
```

So the derivative problem is specifically S-prime, not a broad failure of the empirical S/D sector decomposition.

### 5. Root movement remains tied to the measured residual

Threshold-rank reconstructions give

```text
-DeltaRoot * mean(M') / DeltaM ~= 1
```

on the tested systems, and P45 passes an angular-normalized root-amplitude test. The finite matching root therefore moves through the expected local residual/slope mechanism.

### 6. The derivative sector exposes real subleading structure

Clean 100M full curves resolve that the bare center-slope multiplier is not exactly `2^(3/8)` at current N. On N=185/265:

```text
pure P4[S'] ~ N^-5/4: chi2 = 52.71634 / 2
rank-2/Jordan log:      chi2 =  1.20360 / 2
analytic 1/N:           chi2 =  0.86221 / 2
zero:                    chi2 = 1278.55524 / 2
```

Both predeclared corrections survive. Descriptive chi-square does not reorder frozen chronology.

## Exact finite layer strengthened by the organization pass

Several questions that were previously only proposed now have canonical exact finite results. These do not settle the continuum mechanism, but they improve the foundation from which that mechanism should be derived.

### Typed observable/channel semantics

`cross/either/both/direction`, even/odd combination, probability coordinate, orientation order, raw/normalized status, and scalar/contrast quantity are now explicit descriptors. A scorer must compare identical descriptors or apply a registered exact map.

Both #57 score families are ready under this contract without changing their 2026-08-28 numerical freezes:

- fixed-p primary: exact `either/odd -> cross/odd` identity `D_either=D_cross`;
- intrinsic full-curve cocycle: all primitive observables are `cross`, with size-local P4 normalization before any cross-size comparison.

Broader direction/both semantic audit remains under #146.

### Finite Russo/pivotal identity

For every declared monotone wrapping channel,

```text
M'(p) = pivotal_mass_primal(p) + pivotal_mass_matching(1-p).
```

This is exact on the finite product measure and is independently checked against exact polynomial derivatives on tiny axis, diamond and Gaussian tori; cross also agrees with threshold ranks. The remaining #100 research question is the spin-resolved/continuum pivotal or four-arm structure.

### Exact self-matching tangent

For

```text
p_even = 1/2 + t + lambda
p_odd  = 1/2 + t - lambda,
```

occupation complement acts exactly as `(t,lambda)->(-t,-lambda)`. The N=10 response matrix is

```text
[[0,    0],
 [15/8, 5/4]]
```

for rows `(Rplus,Rminus)` and columns `(t,lambda)`, and independent Bernoulli likelihood-score derivatives reproduce these exact rationals. This gives #61/#155 an explicit microscopic odd tangent; large-N RG decomposition is still open.

### N=26 exact self-matching Beta-family failure

Issue #115 is complete. Exhaustive `2^26` enumeration falsifies both pre-target exact laws with no generalized-Beta rescue fit:

```text
Beta(5,5): first k=5 difference = -96
Beta(7,7): first k=5 difference = +156
```

The new reliability/activation-signature layer explains why N=10 is special. For a centered complement-odd polynomial

```text
M(p)=(2p-1)H(p(1-p)),
```

minimum support plus the minimal-degree bound forces the central-binomial prefix and the Beta law. N=10 lies in that minimal-degree case; N=26 has extra degrees of freedom and is the first nontrivial deformation in this control family.

Threshold ranks can therefore be treated as finite activation/reliability signatures in their own right, linking exact controls, full-curve reconstruction, pivotal mass and profile invariants without invoking a continuum operator.

## What is still weak

### H4 versus H12/H20

Norm-2 `1+i` rotation identifies the odd square-harmonic class, not H4 uniquely. Norm-5 N=325/425 predicts different magnitudes and opposite raw signs for H4 versus H12. This is still the highest-value expensive test.

### Unique radial/operator identity

The `13/8` law is difficult to dismiss over the tested range, but the unique asymptotic exponent and the specific thermal-Q4 assignment remain open. New multipliers, a third full curve, independent controls and modulus fingerprints are more informative than another free fit on the same sizes.

### q=2 versus Jordan correction

N=185/265 rules out the pure S-prime law but not the two frozen corrections. The norm-5 full-curve cocycle, intrinsic thermal coordinates, multi-u response and joint-channel constraints should be used before adding more powers.

### Continuum meaning of matching parity

Empirical pair-exchange sectors are stronger than any proved continuum matching/OPE automorphism. The exact self-matching microscopic tangent narrows the structural question but does not answer how the tangent decomposes into continuum eigenfields.

## Highest-value next work

### 1. Norm-5 Gaussian spectroscopy — #57

The N=325/425 geometry, frozen predictions and both typed scorer entrypoints are ready. Score raw H4 first, then H12/H8/zero in the registered order; reuse the same full curves for the q=2/Jordan functional cocycle. Preserve every child dataset regardless of outcome.

### 2. Third full-curve lineage — #50

Score 145->290 in the frozen order:

1. central residual transfer;
2. bare slope baseline;
3. already-frozen finite-size slope correction;
4. induced root prediction;
5. derivative channels.

### 3. Extract more from existing data

Highest-value low/zero-compute analyses:

- #95 prequential evidence ledger;
- #101 intrinsic quantile-center spectroscopy;
- #119 multi-u functional response;
- #125 joint operator-mixing matrix;
- #122 standardized full threshold distribution;
- #113 local exact-zero analysis;
- #118 dimensionless amplitude-ratio derivation.

For #100 the scalar exact identity is complete; remaining work is orientation-resolved/continuum.

### 4. Choose later expensive geometry by information gain — #102

Do not default to larger N. New geometry/multiplier/modulus leverage should be selected by expected model separation per compute cost.

## Parallel theory/control routes

Keep these independent rather than forcing one story:

- FK/Potts torus-sector derivation (#114);
- pivotal/four-arm H4 geometry (#100/#121);
- torus-modulus and Gaussian-CM/isogeny spectroscopy (#103/#138/#145);
- exactly-critical anisotropy controls (#106);
- configuration-level Euler/Betti lift (#111);
- universal amplitude ratios (#118);
- bounded Galois/reliability/surface-polynomial exact work (#84/#104/#144);
- self-matching tangent RG decomposition (#155).

Pell/modulus production and N=1105 remain gated behind cheaper evidence and an explicit typed/projected lattice observable.

## What to stop doing

Do not spend major effort on:

- interpreting a wrapping-channel sign without a descriptor/exact map;
- adding another correction exponent after every residual;
- generalized Beta fitting after the exact N=26 failure;
- broad PSLQ searches;
- N=1105 before norm-5/new-geometry leverage;
- transfer-matrix width solely for more decimal digits;
- GPU runs without a frozen information target;
- treating multiple derived observables from one histogram as independent evidence.

## Current project thesis

> Square-site/matching finite-size corrections contain reproducible matching-odd and matching-even orientation sectors once exact wrapping-channel semantics are aligned. The central matching-odd residual is compatible with an odd-square-harmonic `N^-13/8` law across independent seeds, exact Gaussian transformations and prospective new geometries, and it moves the finite root through the expected residual/slope mechanism. The matching-even central amplitude also survives the N=185/265 block after exact channel conversion. The exact finite layer now additionally identifies threshold ranks with reliability signatures, the matching slope with pivotal mass, and an explicit self-matching microscopic odd tangent. The remaining hard questions are whether the leading harmonic is genuinely H4, whether the `x=21/4` thermal-family interpretation is correct, and which subleading mechanism explains the prospectively failed pure S-prime law and finite-size slope drift.

That is the line to optimize around until norm-5 or the third full-curve lineage breaks it.
