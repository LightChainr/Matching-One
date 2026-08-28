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
| Central matching-odd `DeltaM` is compatible with `DeltaCos4*N^-13/8` | C3 | P31/P32/P37/P50/P43 | Survives old holdouts, three Gaussian lineages, and prospective N=185/265 geometries |
| N=185/265 x=21/4 H4-like `DeltaM` beats zero and x=17/4 | C3 | P43 / PR #108 | `3.046/2` versus `29.409/2` and `30.246/2` |
| Intrinsic-center `P4[S] ~ N^-1` transfers to N=185/265 | C3 | P48 + prospective score | `1.139/2`, zero `112.540/2` |
| Intrinsic-center `P4[D] ~ N^-13/8` transfers to N=185/265 | C3 | P48 + prospective score | `0.281/2`, zero `29.408/2` |
| Intrinsic-center `P4[D'] ~ N^-5/8` transfers to N=185/265 | C3 | P48 + prospective score | `0.088/2`, zero `59.393/2` |
| Pure `P4[S'] ~ N^-5/4` is sufficient | C3 negative | P48/P43 | Prospectively falsified: `52.716/2` |
| Corrected nonzero `P4[S']` channel exists | C3 | P43/#72 | Both predeclared log and analytic `1/N` corrections survive |
| Issue #43 registered positive DeltaS score is a valid test of cross/even | protocol negative | P43 audit | No: source was `either/even`, target was `cross/even`; original `240.247/2` is a channel-contract failure |
| Exact cross/either map repairs the Issue #43 source channel without target fit | diagnostic | P31/P43 audit | `DeltaS_cross=-DeltaS_either`; repaired diagnostic `0.570/2`; not a retroactive preregistered pass |
| Leading angular harmonic is uniquely H4 | C2 | current data | Not established; #57 norm-5 is the main discriminator |
| Local residual-to-root conversion satisfies `-DeltaRoot*mean(M')/DeltaM ~= 1` | C3 | P35/P45 | Root movement is locally explained by residual and slope |
| Bare finite-size center-slope ratio equals exactly `2^(3/8)` | C2 negative | P49 | False at current precision; finite-size correction resolved |
| Exact C4 self-matching symmetry forces the whole finite polynomial to vanish | C5 negative control | #82 | False; central antisymmetry is exact but polynomial is not identically zero |

## Key prospective numbers

### Matching-odd central sector

```text
N=185 DeltaM = +1.36830e-4 +/- 3.42352e-5
N=265 DeltaM = +1.27110e-4 +/- 3.46783e-5

x=21/4 H4: chi2 = 3.04598 / 2
zero:       chi2 = 29.40938 / 2
x=17/4:     chi2 = 30.24613 / 2
```

This is the strongest new-geometry evidence for the central odd-sector radial law.

### Intrinsic-center P48 four-channel score

Frozen source amplitudes use only N=65,85,130. N=185/265 are independent target geometries and no target amplitude is fit.

```text
P4[S]   ~ N^-1:     chi2 =  1.13878 / 2
P4[D]   ~ N^-13/8:  chi2 =  0.28085 / 2
P4[D']  ~ N^-5/8:   chi2 =  0.08761 / 2
P4[S']  ~ N^-5/4:   chi2 = 52.71634 / 2
```

Only `S'` is a clear pure-law failure.

### S-prime correction score

```text
rank-2/log:    chi2 = 1.20360 / 2
analytic 1/N:  chi2 = 0.86221 / 2
zero:          chi2 = 1278.55524 / 2
```

Both fixed correction models remain viable; descriptive chi-square does not select the mechanism or reorder the frozen chronology.

### Issue #43 registered DeltaS contract audit

The frozen artifact sourced P31 `either/even`, while target production/scoring uses rank-2 `cross/even`.

```text
original registered positive score: chi2 = 240.24721 / 2
exact source-channel map: DeltaS_cross = -DeltaS_either
protocol-repair diagnostic: chi2 = 0.57003 / 2
```

The original score remains a failed preregistration artifact. The repaired value is a lower-status post-reveal protocol diagnostic with zero target refits. It prevents the original failure from being misinterpreted as physical falsification of cross/even N^-1 scaling.

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

The evidence now supports a more coherent picture than the first reading of P43 suggested:

1. **robust matching-odd central sector** — strong finite-size/prospective support for an H4-like `N^-13/8` law;
2. **intrinsic-center parity structure largely survives** — `S`, `D`, and `D'` pure laws transfer to new geometries;
3. **specific derivative correction problem** — `S'` is nonzero but its pure `N^-5/4` law fails and requires finite-size correction;
4. **fixed-coordinate channel bookkeeping matters** — the registered P43 even score mixed an `either/even` source with a `cross/even` target; preserve the failed artifact but do not promote its `240/2` score into a physical sector falsification;
5. **operator theory remains conditional** — unique H4 and a local matching/OPE automorphism are not established.

## Open interpretations

| Interpretation | Level | What would move it forward |
|---|---:|---|
| H4 is the unique leading odd harmonic | C2 | Norm-5 H4-vs-H12 #57 |
| `13/8` is the unique asymptotic radial exponent | C2/C3 | Norm-5 transfer and additional full-curve lineages |
| Minimal relative-`1/N` slope correction explains P49 drift | C0/C2 | Frozen 145->290 full-curve score #50 |
| `S'` correction is ordinary q=2 versus logarithmic/Jordan | C0/C2 | #50/#57 multiplier/full-curve reuse |
| Fixed-coordinate cross/either even observables transport cleanly to intrinsic-center P48 `S` | C0 | #48 thermal-coordinate/channel audit using existing curves |
| `x=21/4` thermal-family spin-4 field is the continuum mechanism | C0/C2 | Unique harmonic, radial competitors, parity controls, corrected derivative structure |
| Matching/complement extends to a local RG/OPE parity automorphism | C0 | #61 plus exact controls |
| `V_<1,4>` explains historical post-`L^-7` behavior | C0 | Axis-annihilator q=3 discrimination |
| Simple PSLQ/algebraic form gives `p_c` | C0 | Low-priority bounded search after provenance constraints |

## Execution priorities

1. **#57 norm-5 N=325/425** — H4-vs-H12 plus radial/correction multiplier leverage. Pilot first, production size from measured power.
2. **#50 N=145->290 full curve** — test the frozen slope correction and induced root prediction on a third lineage.
3. **#48 zero-extra-compute S-prime / coordinate audit** — distinguish q=2 versus log/Jordan and relate fixed-coordinate cross/either observables to intrinsic-center projectors before requesting another run.
4. Tier-B axis-annihilator/exact-zero/control work may continue only when it does not displace the above.

More N=185/265 replicas are not a priority. N=1105 remains gated behind #57/#50.

## Engineering status

- Huawei archive and production tools: canonical on `main`.
- P33 covariance diagnostics: canonical; historical cross-size correlations modest.
- N=185/265 and N=325/425 threshold-rank production: supported.
- Issue #43 primary/secondary scorers and channel-map audit: canonical or in active integration.
- P48 pure-law/correction scoring chain: canonical or in active integration.
- Exact N=1105 projector/minimality and C4 self-matching control: canonical.
- Axis-annihilator Tier-B engine/scorer: canonical.
- CI covers Python 3.9/3.11/3.13 and C++17.

## Explicit non-claims

The project does **not** currently claim:

- a closed form or new exact value for square-site `p_c`;
- proof that `13/8` is the unique asymptotic exponent;
- proof that H4 is the unique harmonic;
- that the original Issue #43 positive DeltaS preregistration passed;
- that its `240/2` failure physically falsifies the matching-even cross sector;
- a successful pure `P4[S'] ~ N^-5/4` law;
- an exact bare `2^(3/8)` finite-size slope ratio;
- a unique q=2 versus Jordan correction mechanism;
- proof of the `x=21/4` LCFT operator identification;
- proof of a full local matching/OPE automorphism.
