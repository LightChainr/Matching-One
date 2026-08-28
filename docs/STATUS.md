# Project Status and Claim Ledger

**Status date:** 2026-08-28

`main` contains the numerical research archive, production tools, exact controls, literature provenance, and active preregistrations. Claim strength is determined by evidence, not by branch location. The execution-facing view is `notes/SYNTHESIS-20260828.md`.

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
| The central matching-odd `DeltaM` sector is compatible with `DeltaCos4*N^-13/8` | C3 | P31/P32/P37/P50/P43 | Survives old holdouts, three Gaussian lineages, and prospective new N=185/265 geometries |
| On N=185/265 the x=21/4 H4-like `DeltaM` law beats zero and the x=17/4 adversary | C3 | P43 / PR #108 | `3.046/2` versus `29.409/2` and `30.246/2` |
| Leading angular harmonic is uniquely H4 | C2 | current data | Not established: `1+i` cannot separate H4 from H12/H20; #57 is the main discriminator |
| Local residual-to-root conversion satisfies `-DeltaRoot*mean(M')/DeltaM ~= 1` | C3 | P35/P45 | Root movement is locally explained by the measured residual and slope |
| Bare finite-size center-slope ratio equals exactly `2^(3/8)` at current sizes | C2 negative | P49 | False at 100M precision; a small finite-size correction is resolved |
| Frozen matching-even `DeltaS ~ +N^-1` assignment is correct | C3 negative | P43 | Prospectively falsified: N=185 and N=265 are both strongly negative |
| Pure `P4[S'] ~ N^-5/4` is sufficient | C3 negative | P48/P43 | Prospectively falsified on new geometries |
| Corrected nonzero `P4[S']` channel exists | C3 | P43/#72 | Both predeclared log and analytic `1/N` corrections survive; mechanism not uniquely selected |
| Exact C4 self-matching symmetry forces the whole finite matching polynomial to vanish | C5 negative control | #82 | False: the central antisymmetry is exact, but the polynomial is not identically zero |

## Key prospective numbers

### N=185/265 matching-odd central sector

```text
N=185 DeltaM = +1.36830e-4 +/- 3.42352e-5
N=265 DeltaM = +1.27110e-4 +/- 3.46783e-5

x=21/4 H4: chi2 = 3.04598 / 2
zero:       chi2 = 29.40938 / 2
x=17/4:     chi2 = 30.24613 / 2
```

This is the strongest new-geometry evidence for the current central odd-sector radial law.

### N=185/265 matching-even central sector

```text
N=185 DeltaS = -6.08154e-5 +/- 8.08957e-6
N=265 DeltaS = -7.02495e-5 +/- 9.38562e-6

frozen positive N^-1: chi2 = 240.24721 / 2
zero:                  chi2 = 112.53891 / 2
```

The simple two-sector conjunction is therefore rejected even though its odd component survives.

### N=185/265 P4[S'] correction score

```text
pure N^-5/4:  chi2 = 52.71634 / 2
rank-2/log:    chi2 =  1.20360 / 2
analytic 1/N:  chi2 =  0.86221 / 2
zero:          chi2 = 1278.55524 / 2
```

Both fixed correction models remain viable. The smaller descriptive statistic is not a license to reorder preregistered chronology.

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

The project should no longer describe the data as one successful pure two-spin-4 model.

A better evidence-respecting decomposition is:

1. **robust matching-odd central sector** — strong finite-size/prospective support for an H4-like `N^-13/8` law;
2. **unresolved matching-even central sector** — the original positive `N^-1` assignment is falsified in sign;
3. **corrected derivative sector** — `S'` is nonzero, but a pure `N^-5/4` law is insufficient and at least one finite-size correction is required;
4. **operator theory remains conditional** — matching/OPE parity and a unique LCFT identification are not established.

## Open interpretations

| Interpretation | Level | What would move it forward |
|---|---:|---|
| H4 is the unique leading odd harmonic | C2 | Norm-5 H4-vs-H12 test #57 |
| `13/8` is the unique asymptotic radial exponent | C2/C3 | Norm-5 transfer and additional full-curve lineages; current new-geometry evidence strongly supports it over the tested range |
| Minimal relative-`1/N` slope correction explains the P49 drift | C0/C2 | Frozen 145->290 full-curve score #50 |
| Matching-even `S` sign reversal is a preasymptotic crossing rather than a different operator sector | C0 | Existing-data signed-sequence analysis #48, then reuse #50/#57 |
| `S'` correction is ordinary q=2 versus logarithmic/Jordan | C0/C2 | Reuse #50/#57 multiplier/full-curve data; no new N=185/265 run |
| `x=21/4` thermal-family spin-4 field is the continuum mechanism | C0/C2 | Unique harmonic, radial competitors, parity controls, and corrected derivative structure |
| Matching/complement extends to a local RG/OPE parity automorphism | C0 | #61 plus exact controls; empirical S/D parity is a weaker statement |
| `V_<1,4>` explains historical post-`L^-7` behavior | C0 | Axis-annihilator experiment; keep conditional until empirical q=3 discrimination |
| Simple PSLQ/algebraic form gives `p_c` | C0 | Low-priority bounded search after provenance constraints |

## Execution priorities

When compute/attention is scarce:

1. **#57 norm-5 N=325/425** — highest-information H4-vs-H12 and radial/correction multiplier test. Start with a small threshold-rank variance pilot, then set production size from measured power.
2. **#50 N=145->290 full curve** — test the already-frozen slope correction and induced root prediction on a third lineage.
3. **#48 zero-extra-compute even/derivative analysis** — explain the `DeltaS` sign reversal and corrected `S'` sequence using existing N=65..265 curves before asking for another run.
4. **Tier-B axis-annihilator/exact-zero/control work** may continue in parallel only when it does not displace the above.

N=1105 remains gated behind #57/#50. More N=185/265 replicas are not a priority; the 500M prospective run is already decisive for the present questions.

## Engineering status

- Huawei archive and production tools: canonical on `main`.
- P33 covariance diagnostics: canonical; historical cross-size correlations are modest and not a global blocker.
- N=185/265 and N=325/425 threshold-rank production: supported on `main`.
- Issue #43 primary/secondary scorers and P48 frozen correction chain: canonical.
- Exact N=1105 projector/minimality and C4 self-matching control: canonical.
- CI covers Python 3.9/3.11/3.13 and C++17.

## Explicit non-claims

The project does **not** currently claim:

- a closed form or new exact value for square-site `p_c`;
- proof that `13/8` is the unique asymptotic exponent;
- proof that H4 is the unique harmonic;
- a successful pure matching-even `N^-1` law;
- a successful pure `P4[S'] ~ N^-5/4` law;
- an exact bare `2^(3/8)` finite-size slope ratio;
- a unique q=2 versus Jordan correction mechanism;
- proof of the `x=21/4` LCFT operator identification;
- proof of a full local matching/OPE automorphism.
