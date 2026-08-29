# Project Status and Claim Ledger

**Status date:** 2026-08-29

`main` is the shared research line. Claim strength is determined by evidence and chronology, not merge status. `docs/RESEARCH-MAP.md` gives relationships among tracks; `notes/SYNTHESIS-20260828.md` gives the current execution view.

## Exact/background facts

| Statement | Level | Status |
|---|---:|---|
| Square-site `p_c` has no known closed form | background | Current project/literature position |
| Square-site and NN+NNN matching-site thresholds satisfy `p_c+p_c_hat=1` | C5 | Exact structural constraint |
| Square-bond and triangular-site thresholds are `1/2` | C5 | Exact controls |
| Threshold-rank histograms are finite activation/reliability signatures under the frozen rank convention | C5 finite | Canonical finite interpretation |
| `M'(p)` equals primal pivotal mass at `p` plus matching pivotal mass at `1-p` | C5 finite | Exact Russo/chain-rule identity |
| The C4 self-matching family has complement tangent `(t,lambda)->(-t,-lambda)` | C5 finite | Exact microscopic `J=-I` |
| Square-bond geometric dual transport swaps primal/dual wrapping and forces `E[D]=0` at `p=1/2` | C5 finite | Exact L=2/L=3 oracle and general finite symmetry argument |

## Strongest current finite-size evidence

| Claim | Level | Evidence | Current interpretation |
|---|---:|---|---|
| Primitive same-`N` Gaussian tori have a nonzero orientation-dependent matching signal | C3 | P31 | Independent-seed confirmation at five frozen sizes |
| The central matching-odd sector is compatible with `DeltaCos4*N^-13/8` | C3 | P31/P32/P37/P50/P43/P57 | Survives held-out/new-geometry and Gaussian multiplier tests |
| N=185/265 favors the frozen x=21/4 H4-like radial law over zero and x=17/4 | C3 | P43 | `3.046/2` vs `29.409/2` and `30.246/2` |
| Norm-5 transfer favors H4 over the frozen H12/H8 aliases | C3 | P57 | H4 `0.4163/2`; H12 `35.1931/2`; H8 `16.0120/2` |
| Norm-5 child block alone proves a nonzero child effect | C3 negative refinement | P57 | No: zero gives `1.77635/2`; the result is a transfer/harmonic discriminator, not a standalone child detection |
| Leading harmonic is globally proven unique H4 | C2 | current data | Not established; tested frozen H12/H8 aliases are strongly disfavored |
| Local root movement obeys `-DeltaRoot*mean(M')/DeltaM ~= 1` | C3 | P35/P45 | Measured residual/slope explains local root motion |
| Frozen matching-even `N^-1` amplitude survives N=185/265 after exact channel conversion | C3 protocol-corrected | P43/#134 | Corrected score `0.5700315436/2` |
| Intrinsic-center `P4[S]`, `P4[D]`, `P4[D']` pure laws survive N=185/265 | C3 | P48 | `1.13878/2`, `0.28085/2`, `0.08761/2` |
| Pure `P4[S'] ~ N^-5/4` is sufficient | C3 negative | P48/P43 | Prospectively falsified: `52.71634/2` |
| q=2 versus rank-2/Jordan is resolved | C2 | P48/P57 | No; norm-5 full-curve score remains inconclusive (`10.648/6` vs `9.020/6`) |
| Simple rank-gap `E[G]=A N^(5/8)+B` correction predicts norm-5 | C2 negative | P57 derived view | Fails at N325/N425: joint `chi2=155.22/2`; source constant-correction fit was already poor |
| A single measured scalar width `w_can` collapses the whole higher thermal jet | C2 negative exploratory | P57 fast jet | Point-level diagnostic fails strongly; do not promote scalar width as the common correction mechanism |

## Key prospective blocks

### P43 — N=185/265 new geometry

```text
N=185 DeltaM = +1.36830e-4 +/- 3.42352e-5
N=265 DeltaM = +1.27110e-4 +/- 3.46783e-5

x=21/4 H4-like: chi2 = 3.04598 / 2
zero:            chi2 = 29.40938 / 2
x=17/4:          chi2 = 30.24613 / 2
```

### P57 — norm-5 N=325/425

500M paired permutations per child; frozen model order H4/H12/H8/zero:

```text
H4:  chi2 =  0.416303764 / 2   p = 0.81208
H12: chi2 = 35.193078878 / 2   p = 2.28e-8
H8:  chi2 = 16.0120      / 2
zero:chi2 =  1.77635     / 2
```

The durable conclusion is narrow but important: the split-prime norm-5 transfer resolves the old `1+i` odd-harmonic alias in favor of the frozen H4 transfer over H12/H8. Because zero is also compatible with the two noisy child measurements, P57 by itself is not a new standalone nonzero-effect detection.

The same raw block gives an inconclusive q=2/Jordan functional score:

```text
q2 analytic:  chi2 = 10.64816 / 6
Jordan/log:   chi2 =  9.02006 / 6
```

These are correlated derived views, not additive evidence beyond the P57 primary block.

## Exact/control and analysis-coordinate progress

### Typed observable semantics

Claim-bearing scores compare identical descriptors or apply a named exact map. The historical Issue #43 even-sector sign error is fixed by

```text
DeltaS_cross = -DeltaS_either.
```

This rule is intentionally narrow: it prevents silent semantic mistakes without blocking exploratory analysis.

### Finite Russo/pivotal identity

```text
M'(p) = pivotal_mass_primal(p) + pivotal_mass_matching(1-p).
```

Exact tiny-system regression is canonical. The open question is the orientation-resolved/continuum pivotal or four-arm structure.

### Exact self-matching and square-bond controls

The C4 self-matching N=10 tangent has response matrix

```text
[[0,    0],
 [15/8, 5/4]].
```

Square-bond geometric dual transport is also canonical: naive bit-complement is not the duality map, while the true transport swaps primal/dual wrapping and forces the odd combination to vanish exactly at `p=1/2`.

### N=26 Beta-family falsification

Exhaustive `2^26` enumeration rejects both frozen laws without a rescue fit:

```text
Beta(5,5): first k=5 coefficient difference = -96
Beta(7,7): first k=5 coefficient difference = +156
```

### Prequential evidence accounting

Issue #95 is complete. One additive primary view is retained per raw random block; roots, slopes, derivatives, score modes, rank gaps and other correlated views may all be analyzed but are not counted as independent votes.

### Krawtchouk/Hermite and rank-gap coordinates

Threshold-rank histograms now support orthonormal Krawtchouk score modes and an exact finite-N Hermite/Krawtchouk generating identity. The joint rank gap

```text
G = K_plus - K_minus
```

has the exact neutral-window area relation

```text
integral_0^1 E[U_{K~Bin(N,p)}] dp = E[G]/(N+1).
```

The simple constant boundary correction to `E[G]` fails on P57. A point-estimate post-reveal test also shows that using the observed scalar `w_can` does not collapse the full higher-mode jet. This shifts mechanism work toward low-rank transfer/mixing rather than another scalar width law.

### Metric-free ratios

The post-P57 definitions

```text
R_I = P4[D'] / (P4[S] * Mbar')
R_T = P4[S'] / (P4[D] * Mbar')
```

cancel the thermal metric. They are useful analysis coordinates, not retrospective P57 predictions. A future held-out score can be frozen against the unrevealed N=145->290 full curve.

## Current interpretation

The evidence-respecting picture is now:

1. **central odd sector** — reproducible and predictive; norm-5 strongly rejects the tested H12/H8 aliases while remaining compatible with H4 transfer;
2. **central even sector** — survives after exact cross/either conversion;
3. **derivative sector** — S, D and D-prime pure laws survive; S-prime requires correction;
4. **simple scalar correction stories are weakening** — both the constant rank-gap boundary model and scalar-width collapse fail on P57 derived views;
5. **compact multivariate transfer is the live mechanism direction** — low-rank matrix/semigroup, operator mixing, and held-out full-curve transfer now have higher information value than adding another scalar exponent;
6. **operator identification remains conditional** — H4 is strongly favored over the tested aliases, but LCFT/Q4 uniqueness and matching/RG parity are not proved.

## Execution priority

1. **#50 N=145->290 full curve** — now the single highest-information new production block. Score the frozen slope/root correction first, then use the same data as a held-out low-rank transfer, Krawtchouk/Hermite, rank-gap and metric-free-ratio test.
2. **Exploit P57 fully with existing data** — covariance-aware higher-mode/low-rank analyses are useful, but they do not require another production run.
3. **#154 norm-4 dyadic closure** — ready if q=2/Jordan or matrix-semigroup ambiguity remains after N290.
4. **#155 self-matching tangent and #159 modulus/Pell controls** — valid parallel work; prioritize whichever gives the largest independent mechanism separation per CPU.
5. **#74 N=1105** — lower current information per cost, not prohibited.

## Engineering/governance status

Research execution is permissive. Only three hard constraints remain: preserve frozen chronology/result history; do not silently score incompatible observable semantics; and do not add correlated views of one raw block as independent primary evidence. Registry/document conflicts do not block scientific integration.

## Explicit non-claims

The project does **not** currently claim:

- a closed form or new exact value for square-site `p_c`;
- proof that `13/8` is the unique asymptotic exponent;
- proof that H4 is globally the unique harmonic/operator;
- proof of a unique q=2 versus Jordan mechanism;
- proof that a scalar width explains the higher thermal-response tower;
- an exact bare `2^(3/8)` finite-size slope ratio;
- proof of the `x=21/4` LCFT operator identification;
- proof of a full local matching/OPE automorphism;
- a rigorous new percolation bound.
