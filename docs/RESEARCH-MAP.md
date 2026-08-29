# Research Map

**Status:** canonical navigation layer  
**Updated:** 2026-08-29

This document explains how the Matching One research program fits together. It does not replace the claim ledger in [`docs/STATUS.md`](STATUS.md), the execution-facing synthesis in [`notes/SYNTHESIS-20260828.md`](../notes/SYNTHESIS-20260828.md), or immutable evidence under `results/`.

The repository has several legitimate organizational axes:

- **scientific question** — topology, threshold signatures, orientation, derivatives, operator identification, exact algebra;
- **evidence block** — P31, P37, P43, P45, P49, P50, exact controls;
- **execution campaign** — for example `results/server-20260828/`;
- **chronology** — frozen predictions, reveals, negative results, errata, later interpretations.

The purpose of this map is to separate those axes so that a later theoretical revision cannot accidentally erase an earlier numerical or exact result.

## 1. Canonical reading order

1. [`README.md`](../README.md) — short external overview.
2. [`docs/STATUS.md`](STATUS.md) — authoritative current claim ledger.
3. [`docs/RESEARCH-MAP.md`](RESEARCH-MAP.md) — relationships among tracks and evidence.
4. [`notes/SYNTHESIS-20260828.md`](../notes/SYNTHESIS-20260828.md) — execution-facing scientific synthesis.
5. [`docs/ROADMAP.md`](ROADMAP.md) — active discriminators and compute priorities.
6. [`analysis/research_ledger.yaml`](../analysis/research_ledger.yaml) — machine-readable questions, evidence blocks, analysis queue, and compute queue.
7. [`analysis/artifact_registry.yaml`](../analysis/artifact_registry.yaml) — machine-readable status of documents, predictions, results, historical protocols, and frontier work.

Old reports, wave notes, and compute queues remain valuable provenance, but they are not competing current-status documents.

## 2. Evidence architecture

Read the project from harder facts upward:

```text
exact matching / finite topology / channel algebra
                    |
                    v
 finite activation signatures / pivotal identities
                    |
                    v
 validated threshold-rank / covariance machinery
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
 Virasoro/LCFT / FK-Potts / pivotal H4 / RG parity / modulus
```

A higher-layer failure does not erase lower-layer evidence. An H12 result would revise H4 without erasing the odd-square-harmonic semigroup evidence; failure of thermal-Q4 would not erase root closure; a different S-prime correction would not erase the prospective failure of pure `N^-5/4`.

## 3. Long-lived research tracks

### A. Exact matching, topology, and observable semantics

**Question:** What is exactly true before asymptotics or CFT enter?

Canonical assets now include:

- finite matching/complement identities;
- arbitrary integer-period torus topology and homology-rank classification;
- Gaussian quotient arithmetic;
- typed observable descriptors for channel, even/odd combination, probability coordinate, orientation order, raw/normalized status, and scalar/contrast quantity;
- fail-closed exact channel mapping, including `DeltaS_cross=-DeltaS_either` and matching-odd `D_either=D_cross` where registered;
- the exact C4 self-matching two-sublattice tangent with microscopic complement action `(t,lambda)->(-t,-lambda)`.

**Boundary:** exact finite pair exchange does not by itself prove a local continuum RG/OPE automorphism.

**Open continuation:** configuration-level Euler/Betti lift (#111), broader direction/both channel audit (#146), and large-N decomposition of the exact self-matching tangent (#155).

### B. Threshold ranks, reliability signatures, and pivotal mass

**Question:** What reusable finite information should expensive simulation preserve?

`K_minus/K_plus` histograms are not merely an implementation trick. Under the frozen rank convention they reconstruct the matching curve and are exact activation/reliability signatures for the finite monotone event. Per-batch histograms and joint moments support later roots, slopes, intrinsic centers, derivative projectors, profiles, and covariance reanalysis without rerunning Monte Carlo.

The finite matching slope also has an exact probability interpretation:

```text
M'(p) = pivotal_mass_primal(p) + pivotal_mass_matching(1-p).
```

This identity is now a canonical tiny-system oracle, including cross-channel agreement with threshold-rank reconstruction.

**Open continuation:** prequential evidence covariance (#95), intrinsic quantile spectroscopy (#101), multi-u response (#119), standardized threshold profiles (#122), joint threshold-gap/signature analysis, and orientation-resolved pivotal/four-arm structure (#100/#121).

### C. Gaussian orientation and the central matching-odd sector

**Question:** Is the microscopic-orientation correction real and predictive at fixed continuum shape/area?

Signal discovery is complete. Same-`N` Gaussian tori, independent seeds, held-out tests, exact Gaussian multiplier lineages, and N=185/265 prospective geometries support a nonzero matching-odd orientation sector.

The compact empirical law over the tested regime is

```text
DeltaM ~ DeltaCos4 * N^(-13/8).
```

**Not established:** H4 rather than H12/H20; unique asymptotic exponent; unique continuum operator.

**Primary discriminator:** #57 norm-5 N=325/425.

### D. Gaussian semigroup, root movement, and thermal metric

**Question:** How does the central sector transform under exact Gaussian multiplication and move the finite root?

P37/P50 provide no-amplitude-fit lineage tests. P35/P45 show local root movement is explained by

```text
DeltaRoot ~= -DeltaM / mean(M').
```

P49 resolves a small but decisive finite-size correction to the bare center-slope multiplier.

**Open continuation:** the frozen 145->290 full-curve correction score (#50) and semigroup residual/correction spectroscopy.

### E. Matching-even and derivative/parity sectors

**Question:** Do `S=(R_G+R_hat)/2` and `D=(R_G-R_hat)/2` expose linked finite-size sectors?

The canonical prospective N=185/265 four-channel score is:

```text
P4[S]   ~ N^-1        survives
P4[D]   ~ N^-13/8     survives
P4[D']  ~ N^-5/8      survives
P4[S']  ~ N^-5/4      pure law fails
```

S-prime is clearly nonzero; analytic q=2 and rank-2/Jordan-log corrections both remain viable.

Both #57 scorer families are now channel-safe without changing their original numerical freezes:

- raw fixed-p primary validates exact `either/odd -> cross/odd` identity before the H4/H12 score;
- intrinsic full-curve cocycle validates identical cross-channel primitive descriptors and size-local P4 normalization before comparing N, 2N, and 5N.

**Open continuation:** #101/#119/#125 and the norm-5 target itself.

### F. Continuum/operator identification

**Question:** What continuum field or probability object produces the measured sector?

The ordinary Virasoro module contains a non-null `c=0,h=5/8` level-4 quasiprimary with bulk `x=21/4`, spin `±4`; in the stated normalization its torus one-point ratio has an exact `g2/E4` fingerprint.

Several independent bridges should remain distinct:

- the exact self-matching microscopic odd tangent -> large-N RG blocks;
- finite Russo identity -> orientation-resolved pivotal/four-arm H4;
- FK/Potts torus-sector Q/thermal derivatives;
- self-dual/self-matching controls;
- torus-modulus one-point fingerprints.

**Boundary:** none of these finite/exact facts yet proves that the lattice matching residual couples uniquely to thermal Q4.

### G. Gaussian CM/isogeny and modular-shape spectroscopy

**Question:** Are multiplier composition and torus shape revealing a representation law rather than only fitted powers?

Gaussian multiplication is an exact finite torus cover/isogeny. Pure scaling fields, finite mixing blocks, and logarithmic/Jordan blocks predict different composition laws. The exact ordinary-Q4 E4 zero supplies an independent modulus axis.

**Gate:** Pell/modulus production remains behind #57/#50 and an operational typed H4-isolating lattice observable. N=1105 remains behind cheaper harmonic/multiplier leverage.

### H. Exact finite matching algebra and reliability structure

**Question:** What exact finite structure exists independently of the continuum story?

This track now contains:

- exact small matching polynomials and roots;
- bounded complex-zero tests, including the failed global imaginary-RMS forecast at exact L=5;
- exact C4 self-matching N=10 control;
- exhaustive N=26 self-matching test that falsifies both frozen `Beta(5,5)` and `Beta(7,7)` laws;
- the self-dual minimal-degree theorem explaining why N=10 `Beta(3,3)` is algebraically special;
- activation/reliability signature and Gaussian/majority shape baselines;
- Boolean/Banzhaf interpretation of the exact self-matching tangent.

These results constrain finite mechanisms; they do **not** imply transcendence or determine the infinite square-site threshold.

### I. Threshold value, literature, and finite-width extrapolation

**Question:** What can be responsibly said about the numerical square-site threshold itself?

Published values remain method-specific evidence with explicit provenance. Blind finite-width work has shown that very accurate next-width prediction can coexist with systematic infinite-intercept drift; Padé/rational alternatives did not automatically cure that drift.

This is a valuable independent methodology track, but it is not the critical path of the orientation/operator program.

## 4. Core evidence chain

| Block | Role | Durable conclusion |
|---|---|---|
| P31 | independent high-stat same-N confirmation | orientation signal is reproducible |
| P32 | held-out radial challenge | fixed `13/8` predicts better than tested zero/free alternatives over tested holdouts |
| P35 | local closure | direct root movement is in the linear residual/slope regime |
| P37 | two fresh norm-2 lineages | parameter-free Gaussian sign/scale relation survives |
| P45 | angular-normalized root amplitude | independently frozen root amplitude is compatible |
| P49 | clean full-curve lineages | central structure survives; bare slope multiplier has a real finite-size correction |
| P50-A | third fixed-p lineage | central no-fit semigroup relation survives a third genealogy |
| P43 | N=185/265 new geometries | matching-odd law beats zero and frozen x=17/4 adversary |
| #134 | exact channel-map erratum | matching-even sign reversal was a source/target semantic bug, not a physical reversal |
| P48 new-geometry score | prospective four-channel test | S, D, D-prime survive; pure S-prime fails |

Do not add these local scores as independent global evidence without shared source/raw-data covariance. #95 exists to build that prequential layer.

## 5. Exact/control evidence added after the original synthesis

| Asset | Status | What it establishes |
|---|---|---|
| Typed channel mapper | exact protocol | channel-bearing scores fail closed unless descriptors match or an exact map exists |
| Russo/pivotal oracle | exact finite | scalar matching slope equals primal+matching pivotal mass |
| Self-matching tangent | exact finite | explicit microscopic pair-exchange-odd tangent and rational N=10 responses |
| N=26 Beta test | exact negative | both frozen finite Beta-family hypotheses fail |
| Reliability signature layer | exact finite | threshold ranks/finite self-dual controls share an activation-signature language |

These strengthen foundations and controls; they do not upgrade the H4/LCFT claim by themselves.

## 6. Current work by resource type

### Existing-data / low-compute

1. #95 prequential evidence ledger.
2. #101 intrinsic quantile-center spectroscopy.
3. #119 multi-u functional response.
4. #125 joint operator-mixing model.
5. #122 standardized full threshold distribution.
6. #113 local exact-zero analysis.
7. #118 dimensionless amplitude-ratio derivation freeze (`R_I`, `R_T`); not a #57 numerical target. Covariance-aware delete-one reconstruction remains.
8. #100/#121 orientation-resolved pivotal/four-arm continuation after the completed scalar Russo identity.

### New production data

1. **#57 norm-5 N=325/425** — H4 versus H12 and linked radial/derivative leverage.
2. **#50 145->290 full curve** — third-lineage slope/root correction test.
3. **#155 large-N self-matching tangent** — only after a source-frozen thermal-orthogonalization rule.
4. **Norm-4 closure** — only if variance/information per CPU justifies it.
5. **Pell/modulus production** — only after an operational H4 projector/observable bridge.
6. **N=1105** — gated behind cheaper tests.

## 7. Negative results are first-class evidence

Preserve and index, rather than hide:

- underpowered large-Pell fixed-p scans;
- wrapping-only GLS `1x` structural negative;
- overfit orientation-difference control-variate pilot;
- Padé failure to cure finite-width signed drift;
- bare finite-size slope multiplier failure;
- prospective pure S-prime failure;
- exact L5 global complex-zero RMS failure;
- graph-cover CRN variance-gain failure;
- exact N26 `Beta(5,5)` and `Beta(7,7)` failures.

## 8. Document/artifact status

Use these labels consistently:

- **canonical-current** — current project truth/navigation;
- **frozen-prediction** — immutable pre-target hypothesis/score contract;
- **immutable-result** — committed evidence archive;
- **topic-derivation** — focused mathematical/theory argument;
- **historical-snapshot** — accurate earlier project state;
- **frontier-open** — unmerged work, not canonical `main` evidence;
- **superseded-protocol** — retained provenance, no longer the active execution entrypoint.

Old queues, notes and reports remain in place; organization is achieved through indexes before file movement.

## 9. Directory contract

```text
constants/      exact/reference constants and relations
data/           literature datasets and provenance
docs/           canonical project-level navigation, status, roadmap
analysis/       machine-readable research/evidence/artifact indexes
notes/          synthesis, theory, derivations, bounded negative results
scripts/        analysis, scorers, exact checks, research utilities
experiments/    frozen or historical protocols and queues
predictions/    preregistered/frozen prediction artifacts
results/        immutable raw and derived research archives
src/            production C++ engines
tests/          scientific-contract, exact-regression, and smoke tests
```

## 10. Update rule

When new evidence lands:

1. preserve the frozen prediction and old result;
2. add/update the evidence block in `analysis/research_ledger.yaml`;
3. update `docs/STATUS.md` only if the current claim boundary changes;
4. update this map only if relationships among tracks change;
5. update `docs/ROADMAP.md` when the information-optimal next discriminator changes;
6. mark old syntheses/queues as historical through the registry rather than rewriting their history.
