# Project Status and Claim Ledger

**Status date:** 2026-08-29

`main` contains the numerical research archive, production tools, exact controls, literature provenance, frozen predictions, protocol errata, and a canonical research/evidence index. Claim strength is determined by evidence and chronology, not by branch location. Use `docs/RESEARCH-MAP.md` for relationships among research tracks and `notes/SYNTHESIS-20260828.md` for the execution-facing scientific view.

## Exact/background facts

| Statement | Level | Status |
|---|---:|---|
| Square-site `p_c` has no known closed form | background | Current project/literature position |
| Square-site and NN+NNN matching-site thresholds satisfy `p_c+p_c_hat=1` | C5 | Exact structural constraint |
| Square-bond and triangular-site thresholds are `1/2` | C5 | Exact controls |
| A rounded decimal is not a definition of `p_c` | governance | Enforced by the literature provenance layer |
| Threshold-rank histograms are exact activation/reliability signatures for the finite monotone event under the frozen rank convention | C5 finite | Canonical finite combinatorial interpretation |
| For the finite matching function, `M'(p)` equals primal pivotal mass at `p` plus matching pivotal mass at `1-p` | C5 finite | Exact Russo/chain-rule identity with independent tiny-torus regression |
| The C4 self-matching two-sublattice family has exact complement tangent `(t,lambda)->(-t,-lambda)` at the center | C5 finite | Exact microscopic `J=-I`; no continuum parity conclusion implied |

## Strongest current finite-size evidence

| Claim | Level | Evidence | Current interpretation |
|---|---:|---|---|
| Primitive same-`N` Gaussian tori have a nonzero orientation-dependent matching signal | C3 | P31 | Independent-seed confirmation at five frozen sizes |
| The central matching-odd `DeltaM` sector is compatible with `DeltaCos4*N^-13/8` | C3 | P31/P32/P37/P50/P43 | Survives old holdouts, three Gaussian lineages, and prospective N=185/265 geometries |
| On N=185/265 the x=21/4 H4-like `DeltaM` law beats zero and x=17/4 | C3 | P43 | `3.046/2` versus `29.409/2` and `30.246/2` |
| Leading angular harmonic is uniquely H4 | C2 | current data | Not established; norm-5 #57 is the main discriminator |
| Local residual-to-root conversion satisfies `-DeltaRoot*mean(M')/DeltaM ~= 1` | C3 | P35/P45 | Root movement is locally explained by the measured residual and slope |
| Frozen matching-even `N^-1` amplitude is compatible with N=185/265 after exact channel conversion | C3 protocol-corrected | P43 + #134 | Frozen source was `either/even`, target was `cross/even`; corrected score `0.5700315436/2` with no refit |
| Intrinsic-center `P4[S]`, `P4[D]`, and `P4[D']` pure laws survive N=185/265 | C3 | P48 new-geometry score | No-refit scores `1.13878/2`, `0.28085/2`, `0.08761/2` |
| Bare finite-size center-slope ratio equals exactly `2^(3/8)` at current sizes | C2 negative refinement | P49 | False at 100M precision; a small finite-size correction is resolved |
| Pure `P4[S'] ~ N^-5/4` is sufficient | C3 negative | P48/P43 | Prospectively falsified: `52.71634/2` |
| Nonzero corrected `P4[S']` channel exists | C3 | P43 | Both predeclared log and analytic `1/N` corrections survive; mechanism not uniquely selected |

## Key prospective numbers

### N=185/265 matching-odd central sector

```text
N=185 DeltaM = +1.36830e-4 +/- 3.42352e-5
N=265 DeltaM = +1.27110e-4 +/- 3.46783e-5

x=21/4 H4-like: chi2 = 3.04598 / 2
zero:            chi2 = 29.40938 / 2
x=17/4:          chi2 = 30.24613 / 2
```

This remains the strongest genuinely new-geometry support for the current central odd-sector radial law.

### N=185/265 matching-even central sector — corrected channel score

The original #108 score compared a frozen `either/even` source prediction to a rank-2 `cross/even` target. Complementary torus topology gives

```text
DeltaS_cross = -DeltaS_either.
```

The no-refit corrected score is

```text
corrected chi2 = 0.5700315436 / 2
marginal z     = +0.6672, -0.1189
```

The original unconverted score remains immutable historical provenance. The executable channel map is now part of the typed observable-semantics layer rather than an inline sign convention.

### Prospective intrinsic-center four-channel score

Using only the frozen N=65/85/130 source amplitudes and the independent N=185/265 target block:

```text
P4[S]   ~ N^-1:     chi2 =  1.13878 / 2
P4[D]   ~ N^-13/8:  chi2 =  0.28085 / 2
P4[D']  ~ N^-5/8:   chi2 =  0.08761 / 2
P4[S']  ~ N^-5/4:   chi2 = 52.71634 / 2
```

This makes the current derivative question narrow: S, D, and D-prime survive their frozen pure laws on new geometries; S-prime is the unique clear pure-law failure among these four channels.

### P4[S'] correction score

```text
pure N^-5/4:   chi2 = 52.71634 / 2
rank-2/log:    chi2 =  1.20360 / 2
analytic 1/N:  chi2 =  0.86221 / 2
zero:           chi2 = 1278.55524 / 2
```

Both fixed correction mechanisms remain viable. Descriptive chi-square alone does not reorder frozen chronology.

## Exact-control progress

### Finite Russo/pivotal identity

For each declared monotone wrapping channel,

```text
M'(p) = pivotal_mass_primal(p) + pivotal_mass_matching(1-p).
```

This is now an exact finite-volume regression oracle on small axis, diamond, and primitive Gaussian tori, and the cross channel also agrees with threshold-rank reconstruction. It gives a probability-theory meaning to the center slope, but it does **not** by itself identify the observed spin-4 correction with a continuum pivotal/four-arm field.

### C4 self-matching tangent

For

```text
p_even = 1/2 + t + lambda
p_odd  = 1/2 + t - lambda,
```

occupation complement gives exactly `(t,lambda)->(-t,-lambda)`. On the exhaustive N=10 control the response matrix, rows `(Rplus,Rminus)` and columns `(t,lambda)`, is

```text
[[0,    0],
 [15/8, 5/4]].
```

Independent Bernoulli likelihood-score derivatives reproduce these rational responses. This is an exact microscopic matching-odd tangent; the large-N decomposition into thermal plus irrelevant RG directions remains open.

### N=26 self-matching Beta-family falsification

The N=10 `Beta(3,3)` threshold law does **not** extend through either of the two frozen N=26 hypotheses. Exhaustive `2^26` enumeration gives:

```text
Beta(5,5): first coefficient difference at k=5 = -96
Beta(7,7): first coefficient difference at k=5 = +156
```

No generalized-Beta rescue fit was performed. The reliability-signature/minimal-degree theorem explains the N=10 result more sharply: `Beta(s,s)` is forced in the minimal-degree self-dual case, while N=26 has extra polynomial degrees of freedom and is the first nontrivial deformation in this control family.

## Current interpretation

The evidence-respecting picture is a two-parity finite-size structure with one unresolved derivative-correction problem:

1. **matching-odd central sector** — strong prospective support for an odd-square-harmonic law compatible with `N^-13/8`;
2. **matching-even central sector** — prospectively compatible with the frozen `N^-1` amplitude after the exact cross/either map;
3. **derivative sector** — S, D, and D-prime pure laws survive N=185/265; S-prime is nonzero but requires a subleading correction;
4. **operator theory remains conditional** — H4 uniqueness, matching/RG parity, and the LCFT family assignment are not proved.

The exact finite-combinatorial layer is now stronger than it was when the first operator notes were written: threshold ranks have a reliability-signature interpretation, the slope has an exact pivotal interpretation, and a self-matching microscopic tangent is explicit. These are bridges and controls; they do not collapse the continuum question into a solved operator identification.

## Open interpretations

| Interpretation | Level | What would move it forward |
|---|---:|---|
| H4 is the unique leading odd harmonic | C2 | Norm-5 H4-vs-H12 test #57 |
| `13/8` is the unique asymptotic radial exponent | C2/C3 | Norm-5 transfer and additional full-curve lineages |
| Minimal relative-`1/N` slope correction explains the P49 drift | C0/C2 | Frozen 145->290 full-curve score #50 |
| `S'` correction is ordinary q=2 versus logarithmic/Jordan | C0/C2 | Norm-5/full-curve cocycle and later multiplier leverage |
| `x=21/4` thermal-family spin-4 field is the continuum mechanism | C0/C2 | Unique harmonic, radial competitors, pivotal/FK bridge, and parity/modulus controls |
| Matching/complement defines a continuum RG parity block | C0/C2 | Large-N continuation of the exact self-matching tangent and independent controls |
| `V_<1,4>` explains historical post-`L^-7` behavior | C0 | Axis-annihilator experiment; keep conditional until q=3 discrimination |
| Simple PSLQ/algebraic form gives `p_c` | C0 | Low-priority bounded search after provenance constraints |

## Execution priorities

When compute/attention is scarce:

1. **#57 norm-5 N=325/425** — highest-information H4-vs-H12 and radial/correction multiplier test. Both the fixed-p primary scorer and the intrinsic full-curve q=2/Jordan scorer now have typed source/target semantic gates on `main`; the frozen numerical contracts themselves were not changed.
2. **#50 N=145->290 full curve** — score the already-frozen slope correction and induced root prediction on a third lineage.
3. **Existing-data analysis** — #95 evidence ledger, #101 quantile-center, #119 multi-u response, #125 joint channel mixing, #122 full threshold profile. For #100, the exact finite Russo stage is complete; the remaining work is orientation-resolved/continuum interpretation.
4. **#102 information-optimal design** chooses any later expensive Gaussian target.

Exact/theory controls (#103/#106/#111/#114/#118/#121) may proceed in parallel. N=1105 and Pell/modulus production remain gated behind cheaper discriminators and an explicit observable bridge.

## Engineering/governance status

- Canonical navigation: `README.md` -> `docs/STATUS.md` -> `docs/RESEARCH-MAP.md` -> `notes/SYNTHESIS-20260828.md` / `docs/ROADMAP.md`.
- Machine-readable knowledge layer: `analysis/research_ledger.yaml` and `analysis/artifact_registry.yaml`.
- Huawei archive and production tools: canonical on `main`.
- P48 N=185/265 four-channel prospective scorer/result: canonical.
- Typed wrapping-channel descriptor/map layer: canonical; a channel mismatch fails closed unless an exact registered map exists.
- #57 fixed-p and full-curve scorer semantic gates: canonical; broader direction/both audit under #146 remains open.
- Exact Russo/pivotal oracle, N=10 self-matching tangent/score-function oracle, N=26 Beta falsification, and reliability-signature tooling: canonical.
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
- that the finite N=26 self-matching negative result determines the infinite square-site threshold;
- a rigorous new percolation bound.
