# Research Map

**Status:** canonical navigation layer  
**Updated:** 2026-08-29

This document explains how the research program fits together. It is not a replacement for the claim ledger in [`docs/STATUS.md`](STATUS.md), the execution-facing synthesis in [`notes/SYNTHESIS-20260828.md`](../notes/SYNTHESIS-20260828.md), or the immutable result archives under `results/`.

The repository has accumulated several valid organizational axes at once:

- **scientific question** — topology, orientation, derivatives, operator identification, exact algebra;
- **experiment/evidence block** — P31, P37, P43, P45, P49, P50, etc.;
- **execution campaign** — for example `results/server-20260828/`;
- **chronology** — frozen predictions, reveals, errata, and later interpretation.

This map separates those axes rather than forcing one directory or issue numbering scheme to carry all of them.

## 1. Canonical reading order

For the current state of the project, read in this order:

1. [`README.md`](../README.md) — short external overview.
2. [`docs/STATUS.md`](STATUS.md) — authoritative current claim ledger.
3. [`docs/RESEARCH-MAP.md`](RESEARCH-MAP.md) — relationships among research tracks and evidence.
4. [`notes/SYNTHESIS-20260828.md`](../notes/SYNTHESIS-20260828.md) — execution-facing scientific synthesis.
5. [`docs/ROADMAP.md`](ROADMAP.md) — active discriminators and compute priorities.
6. [`analysis/research_ledger.yaml`](../analysis/research_ledger.yaml) — machine-readable questions, evidence blocks, analysis queue, and compute queue.
7. [`analysis/artifact_registry.yaml`](../analysis/artifact_registry.yaml) — machine-readable document/result/prediction/frontier registry.

Topic notes and old execution reports remain important, but they are usually **snapshots or derivations**, not competing global status documents.

## 2. Evidence architecture

The project should be read from harder facts upward:

```text
exact matching / finite topology / Gaussian arithmetic
                    |
                    v
       validated measurement machinery
  threshold ranks / K_minus,K_plus / covariance
                    |
                    v
         finite-size empirical structure
 same-N orientation / semigroup / root closure / S,D
                    |
                    v
          mechanism discrimination
 H4 vs H12 / radial corrections / q=2 vs Jordan
                    |
                    v
         continuum interpretation
 Virasoro/LCFT / pivotal / FK-Potts / RG parity / modulus
```

A failure at a higher layer does not erase lower-layer evidence. For example:

- an H12 outcome would revise the H4 interpretation without erasing the measured odd-square-harmonic semigroup structure;
- failure of the thermal-Q4 assignment would not erase the residual-to-root closure;
- a different S-prime correction mechanism would not erase the prospective failure of the pure `N^-5/4` law.

This separation is the main organizing principle for the repository.

## 3. Long-lived research tracks

### A. Exact matching and torus topology

**Question:** What is exactly true before asymptotics or CFT enter?

Core assets include the finite matching identity, complementary wrapping relations, integer-period torus topology, homology-rank classification, exact tiny-system oracles, and Gaussian quotient arithmetic.

This track supplies the semantic and topological contracts used by every numerical result.

**Current boundary:** exact pair/complement identities do not by themselves prove a local RG/OPE automorphism.

**High-value continuation:** typed wrapping-channel algebra, configuration-level Euler/Betti identities, and exact self-matching controls.

### B. Threshold-rank and statistical representation

**Question:** What information should expensive simulations preserve?

The production philosophy is to retain reusable sufficient statistics, especially per-batch integer `K_minus/K_plus` histograms and joint moments, rather than only final roots or fitted amplitudes. These reconstruct `M(p)`, `M'(p)`, intrinsic centers, roots, quantiles, derivative channels, and many later analyses without rerunning Monte Carlo.

**Current strength:** this is a mature research asset, not merely an implementation detail.

**High-value continuation:** evidence-ledger covariance, full threshold-profile analysis, intrinsic quantile spectroscopy, and joint threshold-gap observables.

### C. Gaussian orientation and the central matching-odd sector

**Question:** Is there a real microscopic-orientation correction at fixed continuum torus shape and area?

The signal-discovery phase is complete. Same-`N` Gaussian tori, independent seeds, held-out tests, exact Gaussian multiplier lineages, and N=185/265 prospective geometries all support a nonzero matching-odd orientation sector.

The compact empirical law over the tested regime is

```text
DeltaM ~ DeltaCos4 * N^(-13/8).
```

**Established more strongly than:** any unique CFT interpretation.

**Not yet established:** H4 rather than H12/H20, or uniqueness of the asymptotic radial exponent.

**Primary discriminator:** Issue #57 norm-5 N=325/425.

### D. Gaussian semigroup, root movement, and thermal metric

**Question:** Does the orientation sector transform predictively under exact Gaussian scale/rotation maps, and how does it move the finite root?

P37/P50 supply no-amplitude-fit Gaussian lineage tests. P35/P45 show that local root movement is explained by the measured residual and the center slope:

```text
DeltaRoot ~= -DeltaM / mean(M').
```

P49 additionally resolves a small but statistically decisive finite-size correction to the bare center-slope multiplier.

**Open problem:** determine whether the frozen slope correction predicts a third full-curve lineage and how its correction spectrum relates to the derivative sectors.

### E. Matching-even and derivative/parity sectors

**Question:** Are `S=(R_G+R_hat)/2` and `D=(R_G-R_hat)/2` exposing distinct finite-size fields?

After the exact `either -> cross` correction, the central matching-even `N^-1` amplitude is compatible with N=185/265. The intrinsic-center P48 picture is therefore not a broad parity failure:

```text
P4[S]   ~ N^-1        survives new geometries
P4[D]   ~ N^-13/8     survives new geometries
P4[D']  ~ N^-5/8      survives new geometries
P4[S']  ~ N^-5/4      pure law fails
```

The active derivative problem is narrower: `S'` is clearly nonzero, but the pure leading law is insufficient. Both a predeclared analytic `1/N` correction and a rank-2/Jordan-log correction remain viable.

**High-value continuation:** intrinsic thermal coordinates, multi-u functional response, linked Gaussian multipliers, and joint channel mixing.

### F. Continuum/operator identification

**Question:** What continuum field or probability object produces the measured sector?

The repository contains an exact ordinary-Virasoro result: for `c=0, h=5/8`, a non-null level-4 quasiprimary exists with bulk `x=21/4`, spin `±4`; its torus one-point ratio has an exact weight-4 `g2/E4` fingerprint in the stated normalization.

This makes the thermal-Q4 interpretation concrete, but the key bridge remains open:

```text
lattice matching/topology observable
        ?
continuum thermal-Q4 / logarithmic / defect sector
```

Parallel bridge programs should remain distinct:

- matching/RG tangent-space parity;
- Russo/pivotal and four-arm geometry;
- FK/Potts torus-sector derivatives;
- exact self-dual/self-matching controls;
- torus-modulus one-point fingerprints.

### G. Gaussian CM/isogeny and modular-shape spectroscopy

**Question:** Are Gaussian multipliers and torus moduli exposing a genuine representation law rather than only convenient scaling ratios?

Gaussian multiplication is an exact finite torus cover/isogeny. Pure fields, finite mixing blocks, and logarithmic/Jordan blocks predict different composition laws. The exact thermal-Q4 `E4` result motivates a second, independent shape axis: approach the hexagonal elliptic point, where the ordinary spin-4 one-point vanishes.

**Gate:** do not promote a Pell/modulus production result until the lattice observable/channel projector is explicit and the cheaper norm-5/third-lineage tests are scored.

### H. Exact finite matching algebra and reliability structure

**Question:** What exact structure is present in finite matching polynomials and threshold signatures?

This is an independent exact track. It includes exact small matching polynomials, irreducibility/Galois work, complex-zero falsification tests, self-matching finite controls, and threshold-rank/reliability-signature interpretations.

Negative finite-algebra results constrain simple finite-cell mechanisms but do **not** imply transcendence or absence of some other exact representation for the infinite threshold.

### I. Threshold value, literature, and finite-width extrapolation

**Question:** What can be responsibly said about the numerical square-site threshold itself?

Published threshold values remain method-specific evidence with explicit provenance. The repository does not define `p_c` by one rounded decimal.

Blind finite-width prediction has shown an important methodological distinction: excellent prediction of the next few widths can coexist with systematic drift in the inferred infinite-width intercept. Padé/rational alternatives did not automatically remove that drift.

This track is valuable but is not the critical path of the orientation/operator program.

## 4. Core evidence chain

The current central research chain is:

| Block | Role | Durable conclusion |
|---|---|---|
| P31 | independent high-stat same-N confirmation | orientation-dependent matching signal is reproducible |
| P32 | held-out radial challenge | fixed `13/8` model predicts better than tested zero/free alternatives over the tested regime |
| P35 | local closure | direct root movement is in the linear residual/slope regime |
| P37 | two fresh norm-2 lineages | parameter-free Gaussian sign/scale relation survives |
| P45 | angular-normalized root amplitude | root amplitude is compatible with an independently frozen prediction |
| P49 | clean full-curve lineages | central structure survives; bare slope multiplier has a real finite-size correction |
| P50-A | third fixed-p lineage | the central no-fit semigroup relation survives a third genealogy |
| P43 | N=185/265 new geometries | matching-odd law beats zero and the frozen x=17/4 adversary |
| #134 | exact channel-map erratum | matching-even sign reversal was a source/target channel mismatch, not a physical reversal |
| P48/P43 S-prime | prospective derivative test | pure `N^-5/4` fails; q=2 and Jordan/log corrections remain live |

Do not add these local scores as independent global evidence without accounting for shared source parameters and raw-data covariance. Issue #95 exists to build that prequential evidence layer.

## 5. Current work by resource type

### Existing-data / low-compute analysis

Highest-value analyses that mostly reuse committed full curves or exact artifacts:

1. **#95 prequential evidence ledger** — chronology-locked predictive scores and evidence coverage.
2. **#101 intrinsic quantile-center spectroscopy** — isolate nonlinear thermal-coordinate effects.
3. **#119 multi-u thermal response** — use the frozen `u={0,.025,.05}` vector as one correlated functional observation.
4. **#125 joint operator-mixing model** — force one physical correction to move several channels coherently.
5. **#100 pivotal/Russo bridge** — exact finite-volume identity first, then orientation-resolved pivotal structure.
6. **#122 standardized full threshold distribution** — use the histogram asset beyond local derivatives.
7. **#113 local exact-zero analysis** — use existing exact roots near the physical root rather than global root-cloud summaries.
8. **#118 universal amplitude-ratio derivation** — derive dimensionless combinations before measuring them.

### New production data

Order by expected discrimination value rather than by issue number:

1. **#57 norm-5 N=325/425** — H4 versus H12 and linked radial/derivative leverage.
2. **#50 145->290 full curve** — third-lineage slope/root correction test.
3. **Self-matching / self-dual tangent controls** — test parity and RG-direction interpretations independently.
4. **Norm-4 closure** — only if variance/power shows it is competitive for q=2 versus Jordan.
5. **Pell/modulus production** — gated on typed channel semantics and an operational H4-isolating observable.
6. **N=1105** — gated behind cheaper harmonic/multiplier information.

## 6. Negative results are first-class evidence

Preserve and index, rather than hide, at least these results:

- large-Pell fixed-p scans were underpowered;
- wrapping-only GLS gave no variance reduction because several differences are configuration-identical;
- pilot-fitted orientation-difference control variates overfit;
- Padé did not cure the signed finite-width extrapolation drift;
- the bare finite-size slope multiplier is not exactly `2^(3/8)` at current sizes;
- the pure `P4[S'] ~ N^-5/4` law fails prospectively;
- the global complex-zero imaginary-RMS forecast failed at exact L=5;
- an exact graph-cover CRN construction did not provide useful variance reduction in its pilot;
- simple finite self-matching Beta-family extensions are exact-control hypotheses, not assumptions to preserve after a frozen failure.

A negative result should retain its original protocol and raw artifact even when a better method or interpretation follows.

## 7. Document and artifact status

Use these labels consistently:

- **canonical-current** — current project-level truth or navigation (`README`, `STATUS`, `RESEARCH-MAP`, `SYNTHESIS`, `ROADMAP`);
- **frozen-prediction** — pre-target hypothesis/score contract; never silently rewrite after reveal;
- **immutable-result** — raw/derived evidence archive and report for a completed run;
- **topic-derivation** — mathematical/theoretical derivation whose local claims may remain valid even if global interpretation changes;
- **historical-snapshot** — synthesis, queue, or report that accurately records an earlier project state;
- **frontier-open** — work completed or proposed on an unmerged branch/PR; do not cite it as canonical `main` evidence;
- **superseded-protocol** — retained for provenance but no longer the execution entry point.

Old server queues, wave notes, and reports are not to be deleted merely because priorities changed. Their role is provenance.

## 8. Directory contract

```text
constants/      exact/reference constants and relations
data/           literature datasets and provenance
docs/           canonical project-level navigation, status, roadmap
analysis/       machine-readable research/evidence/artifact indexes
notes/          synthesis, theory, derivations, bounded negative results
scripts/        analysis, scorers, exact checks, research utilities
experiments/    frozen or historical experiment protocols and queues
predictions/    preregistered/frozen prediction artifacts
results/        immutable raw and derived research archives
src/            production C++ engines
tests/          scientific-contract, exact-regression, and smoke tests
```

The directory contract is intentionally conservative: organization is achieved through indexes and typed relations before any mass file movement.

## 9. Observable semantics

Every channel-bearing comparison should identify, at minimum:

- topology channel (`cross`, `either`, `both`, direction);
- primal/matching or even/odd combination;
- probability coordinate (`p` or complement);
- signed orientation order;
- raw contrast versus angular-normalized projector;
- scalar value versus orientation contrast.

A scorer must compare identical descriptors or apply an exact registered map. Hash equality is not a substitute for semantic equality.

## 10. Update rule

When new evidence lands:

1. preserve the frozen prediction and old result;
2. add the new evidence block to `analysis/research_ledger.yaml`;
3. update `docs/STATUS.md` only if the current claim boundary changes;
4. update this map only if the relationship among tracks changes;
5. update `docs/ROADMAP.md` when the information-optimal next discriminator changes;
6. mark old synthesis/queue documents as historical through the registry rather than rewriting their history.

This keeps the repository cumulative: new interpretation should connect to old evidence, not overwrite it.
