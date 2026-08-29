# Research Synthesis — 2026-08-29 execution sync

This is the execution-facing view of Matching One. It is deliberately shorter and more permissive than earlier wave notes: useful analysis, exact work and pilots may proceed by default. `docs/STATUS.md` controls claim language; `docs/ROADMAP.md` ranks information gain.

## Current thesis

The strongest finite-size picture remains a **two-parity orientation structure with one unresolved derivative/metric mechanism**:

```text
matching-odd central sector
    -> robust and prospectively reproduced
    -> compatible with DeltaCos4 * N^-13/8
    -> root motion locally follows DeltaRoot ~= -DeltaM/M'

matching-even central sector
    -> compatible with the frozen N^-1 amplitude after exact either->cross conversion

intrinsic derivative sectors
    -> S, D and D-prime pure laws survive N=185/265
    -> S-prime is decisively nonzero
    -> pure N^-5/4 fails prospectively
    -> q=2 and Jordan-log corrections both remain live
```

The main uncertainty is no longer whether an orientation signal exists. It is **what compact mechanism transports the full response across size, multiplier and geometry**.

## Durable evidence

### Central matching-odd signal

Independent N=65,85,130,145,170 data and prospective Gaussian lineages are compatible with

```text
DeltaM ~ DeltaCos4 * N^-13/8.
```

The N=185/265 new-geometry block gives

```text
N=185 DeltaM = +1.36830e-4 +/- 3.42352e-5
N=265 DeltaM = +1.27110e-4 +/- 3.46783e-5

x=21/4 H4-like: chi2 = 3.04598 / 2
zero:            chi2 = 29.40938 / 2
x=17/4:          chi2 = 30.24613 / 2
```

The third fixed-coordinate 145->290 lineage also passed prospectively with residual z about `-0.48`.

This is strong C3 finite-size evidence. H4 uniqueness, unique asymptotic exponent and unique continuum operator remain open.

### Matching-even correction

The historical N=185/265 sign conflict was a source/target observable mismatch. Exact torus topology gives

```text
DeltaS_cross = -DeltaS_either.
```

With that conversion and no refit,

```text
corrected chi2 = 0.5700315436 / 2
z = +0.667, -0.119.
```

The old unconverted score remains provenance; it is excluded from current evidence.

### Derivative sector

Using frozen source amplitudes and independent N=185/265 targets:

```text
P4[S]   ~ N^-1:     chi2 =  1.13878 / 2
P4[D]   ~ N^-13/8:  chi2 =  0.28085 / 2
P4[D']  ~ N^-5/8:   chi2 =  0.08761 / 2
P4[S']  ~ N^-5/4:  chi2 = 52.71634 / 2
```

For S-prime:

```text
pure N^-5/4:   chi2 = 52.71634 / 2
rank-2/log:    chi2 =  1.20360 / 2
analytic 1/N:  chi2 =  0.86221 / 2
zero:           chi2 = 1278.55524 / 2
```

So the open mechanism question is narrow: the channel is real, the pure law fails, and two fixed corrections remain viable.

## The finite layer is now an analysis engine

The repository now preserves enough structure that one expensive threshold-rank block can answer many questions without rerunning simulation.

### Exact/semantic primitives

- typed channel/combination/coordinate/order/normalization descriptors;
- exact registered cross/either maps;
- arbitrary integer-period torus backend with exact winding arithmetic;
- exact finite Russo/pivotal identity;
- exact C4 self-matching tangent;
- exact N=26 Beta-family falsification and reliability-signature interpretation.

### Evidence accounting

The prequential ledger from #95 is complete. One raw random block contributes one additive primary view; correlated root/slope/derivative/score-mode diagnostics remain available without becoming duplicate evidence votes.

### Krawtchouk/Hermite response coordinates

The same threshold-rank histograms can be projected onto orthonormal binomial Krawtchouk score modes. Mode 0 and mode 1 exactly recover value and first-derivative information in new coordinates, while modes `r>=2` expose a higher thermal-response vector.

The exact finite-N generating function is

```text
R_N(p) = sum_r c_r sqrt(C(N,r))
         ((p-p0)/sqrt(p0(1-p0)))^r.
```

This means the score tower is one full-curve Taylor jet, not a collection of unrelated observables.

### Paired rank-gap width

The joint statistic

```text
G = K_plus - K_minus
```

contains information not reconstructible from the two marginal rank laws. Exactly,

```text
integral_0^1 E[U_{K~Bin(N,p)}] dp = E[G]/(N+1),
```

where `U` is the neutral-window indicator. This provides a canonical width coordinate

```text
w_can(N) = N^(3/8) E[G]/(N+1).
```

The source fit and its N=325/425 target predictions are frozen before the norm-5 reveal. The strongest immediate mechanism test is therefore width-only collapse of the Hermite/Krawtchouk jet, followed by q=2 versus Jordan cocycle scoring if a residual direction remains.

### Low-rank matrix/semigroup discovery

The retrospective P48/full-curve analysis now has a compact matrix-state candidate: a rank-1 identity block plus a rank-2 thermal block gives a much better organizing picture than assigning a new independent exponent to every scalar channel. Its Jordan-vs-ordinary ranking is post-reveal discovery, not prospective evidence, but it creates a sharper held-out transfer test for N=290 and norm-5.

### Information-oriented geometry design

The maximin Gaussian design tooling from #139 is canonical. It confirms norm-5 as a high-leverage H4/H12/radial discriminator under its planning model and supplies exact multiplier/harmonic arithmetic for choosing later geometries by model separation per compute cost.

## Highest-value execution

### 1. #57 N=325/425 — one block, many questions

Do not treat norm-5 as a single H4/H12 number. Produce full threshold-rank histograms and joint moments, then score the same raw block in this order:

1. frozen H4/H12/H8/zero fixed-p primary;
2. intrinsic/full-curve q=2 versus Jordan cocycle;
3. frozen Krawtchouk thermal vector `(D2,S3,D4,S5,D6)`;
4. frozen paired rank-gap N=325/425 targets;
5. canonical rank-gap width -> Hermite/Krawtchouk width-collapse residual;
6. root/slope closure, derivative, multi-u and low-rank-transfer diagnostics;
7. update the prequential ledger with one additive primary view for the raw block.

A small pilot or exploratory block may start whenever useful. If optional-stopping guarantees are desired, use the frozen sequential rule from the beginning; otherwise the data remain ordinary exploratory/fixed-count evidence.

### 2. #50 N=145->290 full curve

The fixed-coordinate result is already prospective. The full curve should test:

1. central residual transfer;
2. bare center-slope multiplier;
3. frozen finite-size slope correction;
4. direct and induced root transfer;
5. derivative channels;
6. Krawtchouk/rank-gap coordinates;
7. held-out low-rank matrix/semigroup transfer.

This is a third-lineage mechanism test, not another signal-discovery run.

### 3. Existing-data analysis continues in parallel

Run useful work without waiting for a roadmap promotion:

- #101 intrinsic quantile centers;
- #119 multi-u functional response;
- #125 joint operator mixing;
- #122 standardized threshold profile;
- #100/#121 orientation-resolved pivotal/four-arm structure;
- #113 exact local-zero structure;
- #118 dimensionless amplitude ratios;
- #180 low-rank full-curve transfer operator.

## Later compute is a priority choice, not a gate

Norm-4 N=260/340, large-N self-matching tangent, Pell/modulus geometry and N=1105 are all valid exploratory directions. Their present ordering reflects expected information per CPU, not permission.

Use the maximin/information-design tools to decide when one of them outranks extending #57 or #50. A strong operator-specific Pell/modulus claim still needs an explicit observable bridge, but generating exploratory data does not.

## Parallel theory/control routes

Keep multiple bridges alive until one creates a sharper score:

- FK/Potts torus-sector derivation;
- pivotal/four-arm anisotropy;
- torus-modulus and Gaussian-CM/isogeny spectroscopy;
- exactly-critical anisotropy controls;
- Euler/Betti/topological polynomial identities;
- universal amplitude ratios;
- self-matching RG tangent;
- bounded exact/Galois work.

Theory earns priority when it produces an independent derivation, exact control, or parameter-free/low-dimensional target that current data can score.

## Avoid low-information loops

Low priority rather than forbidden:

- another free exponent fit on the old five sizes;
- more N=185/265 replicas merely repeating an already completed test;
- generalized Beta rescue fitting after the exact N=26 failure;
- broad PSLQ without a new structural constraint;
- transfer-matrix width solely for threshold digits;
- large production that stores only a final scalar instead of reusable sufficient statistics;
- treating several views from one histogram as independent evidence.

## Current project thesis

> Square-site/matching finite-size corrections contain reproducible matching-odd and matching-even orientation sectors once exact wrapping-channel semantics are aligned. The central matching-odd residual is compatible with an odd-square-harmonic `N^-13/8` law across independent seeds, exact Gaussian transformations and prospective new geometries, and it moves the finite root through the expected residual/slope mechanism. The matching-even central amplitude survives the N=185/265 block after exact channel conversion. Threshold ranks now also provide a reusable finite analysis state: reliability signatures, pivotal mass, Krawtchouk/Hermite response modes, paired rank-gap width and low-rank transfer structure. The next decisive step is to use N=325/425 and N=290 full curves as multi-question held-out blocks that discriminate harmonic content, thermal metric and compact transfer mechanisms at once.
