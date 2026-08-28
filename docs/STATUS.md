# Project Status and Claim Ledger

**Status date:** 2026-08-28

`main` contains the numerical research archive, production tools, exact controls, literature provenance, preregistered predictions, and protocol errata. Claim strength is determined by evidence and chronology, not by branch location. The execution-facing view is `notes/SYNTHESIS-20260828.md`.

## Exact/background facts

| Statement | Level | Status |
|---|---:|---|
| Square-site `p_c` has no known closed form | background | Current project/literature position |
| Square-site and NN+NNN matching-site thresholds satisfy `p_c+p_c_hat=1` | C5 | Exact structural constraint |
| Square-bond and triangular-site thresholds are `1/2` | C5 | Exact controls |
| A rounded decimal is not a definition of `p_c` | governance | Enforced by the literature provenance layer |

## Strongest current finite-size evidence

| Claim | Level | Evidence | Current interpretation |
|---|---:|---|---|
| Primitive same-`N` Gaussian tori have a nonzero orientation-dependent matching signal | C3 | P31 | Independent-seed confirmation at five frozen sizes |
| The central matching-odd `DeltaM` sector is compatible with `DeltaCos4*N^-13/8` | C3 | P31/P32/P37/P50/P43 | Survives old holdouts, three Gaussian lineages, and prospective N=185/265 geometries |
| On N=185/265 the x=21/4 H4-like `DeltaM` law beats zero and x=17/4 | C3 | P43 / #108 | `3.046/2` versus `29.409/2` and `30.246/2` |
| Frozen matching-even `N^-1` amplitude is compatible with N=185/265 after exact channel conversion | C3 protocol-corrected | #108 + #134 | Frozen source was `either/even`, target was `cross/even`; corrected score `0.5700/2` with no refit |
| Intrinsic-center `P4[S] ~ N^-1` transfers to prospective N=185/265 | C3 | P48 prospective score | `1.1388/2`, zero `112.540/2` |
| Intrinsic-center `P4[D] ~ N^-13/8` transfers to prospective N=185/265 | C3 | P48 prospective score | `0.2809/2`, zero `29.408/2` |
| Intrinsic-center `P4[D'] ~ N^-5/8` transfers to prospective N=185/265 | C3 | P48 prospective score | `0.0876/2`, zero `59.393/2` |
| Pure `P4[S'] ~ N^-5/4` is sufficient | C3 negative | P48/P43 | Prospectively falsified: `52.7163/2` |
| Nonzero corrected `P4[S']` channel exists | C3 | P43/#72 | Both predeclared log and analytic `1/N` corrections survive; mechanism not uniquely selected |
| Leading angular harmonic is uniquely H4 | C2 | current data | Not established; norm-5 #57 is the main discriminator |
| Local residual-to-root conversion satisfies `-DeltaRoot*mean(M')/DeltaM ~= 1` | C3 | P35/P45 | Root movement is locally explained by the measured residual and slope |
| Bare finite-size center-slope ratio equals exactly `2^(3/8)` at current sizes | C2 negative refinement | P49 | False at 100M precision; a small finite-size correction is resolved |
| Exact C4 self-matching symmetry forces the whole finite matching polynomial to vanish | C5 negative control | #82 | False: central antisymmetry is exact, polynomial is not identically zero |

## Key prospective numbers

### N=185/265 matching-odd central sector

```text
N=185 DeltaM = +1.36830e-4 +/- 3.42352e-5
N=265 DeltaM = +1.27110e-4 +/- 3.46783e-5

x=21/4 H4-like: chi2 = 3.04598 / 2
zero:            chi2 = 29.40938 / 2
x=17/4:          chi2 = 30.24613 / 2
```

This is the strongest new-geometry evidence for the current central odd-sector radial law.

### N=185/265 matching-even central sector — corrected channel score

The original #108 score compared a frozen `either/even` source prediction to a rank-2 `cross/even` target. Complementary torus topology gives

```text
DeltaS_cross = -DeltaS_either.
```

The original frozen either means become the cross predictions

```text
N=185 predicted cross DeltaS = -6.75216374588e-5
N=265 predicted cross DeltaS = -6.89194469703e-5
```

against observations

```text
N=185 observed = -6.08153762334e-5 +/- 8.08956556558e-6
N=265 observed = -7.02495078452e-5 +/- 9.38562007761e-6
```

Using the unchanged fully correlated source-amplitude uncertainty:

```text
corrected chi2 = 0.5700315436 / 2
marginal z     = +0.6672, -0.1189
```

This is a protocol correction with zero target-fit parameters. The original #108 score remains historical provenance; PR #134 supplies the erratum and regression.

### Prospective intrinsic-center P48 four-channel score

The source amplitudes are frozen from N=65,85,130. N=185/265 are independent targets and no target parameter is fit.

```text
P4[S]   ~ N^-1:     chi2 =  1.13878 / 2   zero = 112.53974 / 2
P4[D]   ~ N^-13/8:  chi2 =  0.28085 / 2   zero =  29.40844 / 2
P4[D']  ~ N^-5/8:   chi2 =  0.08761 / 2   zero =  59.39319 / 2
P4[S']  ~ N^-5/4:   chi2 = 52.71634 / 2   zero = 1278.55524 / 2
```

The empirical P48 parity pattern is therefore prospectively stable in `S`, `D`, and `D'`. `S'` is the unique clear pure-law failure among these four channels.

### N=185/265 P4[S'] correction score

```text
rank-2/log:    chi2 = 1.20360 / 2
analytic 1/N:  chi2 = 0.86221 / 2
zero:          chi2 = 1278.55524 / 2
```

Both fixed correction models remain viable. Descriptive chi-square alone does not reorder frozen chronology.

### Earlier orientation/root evidence

```text
P31 pooled A4 = 0.7885 +/- 0.0352
P37 65->130   = -0.31382 +/- 0.0908
P37 85->170   = -0.34095 +/- 0.1118
frozen norm-2 ratio = -0.3242098887...
P50 145->290 fixed-p residual z = -0.483
```

P45 angular-normalized root amplitudes:

```text
A_p(65) = 0.42034 +/- 0.02157
A_p(85) = 0.39495 +/- 0.03078
frozen  = 0.45101 +/- 0.02013
```

P49 clean slope result:

```text
observed ratios = 1.2939835, 1.2943776
bare target     = 1.2968396
bare slope chi2 = 6412.89 / 2
```

## Current interpretation

The best evidence-respecting picture is a two-parity finite-size structure with one specifically unresolved derivative channel:

1. **matching-odd central sector** — strong prospective support for an odd-square-harmonic law compatible with `N^-13/8`;
2. **matching-even central sector** — prospectively compatible with the frozen `N^-1` amplitude after exact cross/either channel mapping;
3. **intrinsic-center parity pattern** — `S`, `D`, and `D'` pure laws independently transfer to N=185/265;
4. **S-prime correction problem** — `S'` is nonzero but pure `N^-5/4` fails and requires a subleading correction;
5. **operator theory remains conditional** — H4 uniqueness, matching/OPE parity, and LCFT family assignment are not proved.

The main lesson from #134 is methodological: wrapping-channel semantics are part of the prediction. The main lesson from the prospective P48 score is scientific: the channel erratum does not merely rescue one central amplitude; most of the intrinsic-center parity spectrum also survives independently.

## Open interpretations

| Interpretation | Level | What would move it forward |
|---|---:|---|
| H4 is the unique leading odd harmonic | C2 | Norm-5 H4-vs-H12 test #57 |
| `13/8` is the unique asymptotic radial exponent | C2/C3 | Norm-5 transfer and additional full-curve lineages |
| Minimal relative-`1/N` slope correction explains the P49 drift | C0/C2 | Frozen 145->290 full-curve score #50 |
| `S'` correction is ordinary q=2 versus logarithmic/Jordan | C0/C2 | #50/#57 multiplier/full-curve leverage; multi-u/joint-channel analysis |
| `x=21/4` thermal-family spin-4 field is the continuum mechanism | C0/C2 | Unique harmonic, radial competitors, pivotal/FK bridge, and parity controls |
| Matching/complement extends to a local RG/OPE parity automorphism | C0 | Exact controls + theory; empirical pair-exchange parity is weaker |
| `V_<1,4>` explains historical post-`L^-7` behavior | C0 | Axis-annihilator experiment; keep conditional until q=3 discrimination |
| N=10 self-matching `Beta(3,3)` extends to a finite exact family | C0 | Pre-frozen N=26 exact falsification #115 |
| Simple PSLQ/algebraic form gives `p_c` | C0 | Low-priority bounded search after provenance constraints |

## Execution priorities

When compute/attention is scarce:

1. **#57 norm-5 N=325/425** — highest-information H4-vs-H12 and radial/correction multiplier test.
2. **#50 N=145->290 full curve** — score the already-frozen slope correction and induced root prediction on a third lineage.
3. **Zero/low-compute analysis** — #48 S-prime correction/channel-coordinate audit, #95 evidence ledger, #100 pivotal bridge, #101 quantile-center, #119 multi-u response, #125 joint channel mixing.
4. **#102 information-optimal design** chooses any later expensive Gaussian target.

Exact/theory controls (#103/#106/#111/#114/#115/#118/#121) may proceed in parallel. N=1105 remains gated behind cheaper multipliers and clear information need.

## Engineering/governance status

- Huawei archive and production tools: canonical on `main`.
- N=185/265 and N=325/425 threshold-rank production: supported on `main`.
- Issue #43 primary/secondary scorers and P48 frozen correction chain: canonical.
- Cross/either Issue #43 protocol erratum: canonical via PR #134.
- Prospective P48 four-channel scorer/result: canonical once its current integration PR lands.
- Exact N=1105 projector/minimality and C4 self-matching control: canonical.
- CI covers Python 3.9/3.11/3.13 and C++17.

## Explicit non-claims

The project does **not** currently claim:

- a closed form or new exact value for square-site `p_c`;
- proof that `13/8` is the unique asymptotic exponent;
- proof that H4 is the unique harmonic;
- proof of a unique q=2 versus Jordan correction mechanism;
- an exact bare `2^(3/8)` finite-size slope ratio;
- proof of the `x=21/4` LCFT operator identification;
- proof of a full local matching/OPE automorphism;
- a rigorous new percolation bound.
