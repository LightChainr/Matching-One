# Post-confirmation review: no-fit doubling passes; full-curve/root scaling remains open

Status: review of PR #21 at server head `6d2d68a62b433e337b97eadaf1870cb58e2f7666`.

This note supersedes the checkpoint-specific conclusions in earlier versions. The square-site matching-orientation effect has now passed three distinct finite-size tests: an independent seed, a held-out radial/harmonic challenge, and a parameter-free exact Gaussian doubling relation. The remaining questions are full-curve thermal projection, cross-size covariance, the `N^-2` root law, higher odd harmonics, clean production provenance, variance reduction, and operator identification.

## 1. The same-N orientation effect is confirmed on the tested finite sizes

The frozen P31 replay used an independent seed, disjoint counters, 100 equal batches and 100,000,000 paired replicas at every prescribed size:

| `N` | pair | `Delta M` | batch SE | z | `N^(13/8) Delta M / Delta cos(4 theta)` |
|---:|---|---:|---:|---:|---:|
| 65 | `(8,1)/(7,4)` | `+1.24948e-3` | `7.80e-5` | 16.03 | 0.8093 |
| 85 | `(9,2)/(7,6)` | `+1.01189e-3` | `9.01e-5` | 11.23 | 0.8666 |
| 130 | `(11,3)/(9,7)` | `+4.67000e-4` | `8.95e-5` | 5.22 | 0.9330 |
| 145 | `(12,1)/(9,8)` | `+4.42250e-4` | `8.38e-5` | 5.27 | 0.7501 |
| 170 | `(13,1)/(11,7)` | `+2.37640e-4` | `9.20e-5` | 2.58 | 0.6277 |

Every sign agrees with `Delta cos(4 theta)`. The seed is statistically compatible with the earlier 30-million campaign at overlapping sizes. Pooling seeds within each size gives the descriptive common-amplitude summary

\[
A_4=0.7885\pm0.0352,
\qquad
\chi^2=1.53\quad(4\ \mathrm{dof}).
\]

The existence, sign, and approximate radial collapse of the finite-size orientation effect are no longer exploratory claims.

## 2. The preregistered H4 model passes the first held-out challenge

P32 trained only on `N=65,85,130`; `N=145,170` were hidden until scoring. Every available independent seed row was retained. The held-out scores are:

| model | held-out chi-square |
|---|---:|
| fixed `13/8`, H4 | 1.058 |
| fixed `13/8` plus one power correction | 1.712 |
| fixed `13/8` plus logarithmic amplitude | 1.726 |
| fixed `13/8`, H4+H8 | 1.100 |
| free radial exponent | 1.661 |
| zero effect | 37.32 |

The simplest frozen law

\[
\boxed{
\Delta M_N=A_4\,\Delta\cos(4\theta)\,N^{-13/8}
}
\]

has the best held-out score and conditioning. The fitted even-parity higher harmonic is

\[
A_8=-0.0345\pm0.0542,
\]

so H8 is unresolved and does not improve prediction. The power, logarithmic, and free-exponent alternatives are not required at current precision and are poorly conditioned relative to H4.

This is a successful Stage-1 structural test. It favors H4 with the preregistered `13/8` power, but does not uniquely establish the asymptotic exponent.

## 3. P37 passes an exact, parameter-free doubling test

A stronger test was preregistered in commit `ceb7a578eaface6f6882278c567b4d893014f59c` before the fresh production run. Multiplication of a Gaussian period by `1+i` gives

\[
N\mapsto2N,
\qquad
\theta\mapsto\theta+\pi/4.
\]

For a leading `cos(4 theta) N^-13/8` term, both the angular phase and radial magnitude are fixed:

\[
\boxed{
\frac{\Delta M_{2N}}{\Delta M_N}=-2^{-13/8}
=-0.3242098886627524.
}
\]

No amplitude, exponent, or critical probability is fitted. The fresh seed uses clean source commit `80fbdd1e9a380a87a3c56dec7795ceebb0ada23e`, seed `2026100101`, counters `[2000000000,2100000000)`, 100 equal batches, and 100,000,000 replicas at each of `N=65,85,130,170`.

Lineage order, rather than display order, gives:

```text
65 -> 130:
  (8,1) -> (9,7)
  (7,4) -> (11,3)

85 -> 170:
  (9,2) -> (11,7)
  (7,6) -> (13,1)
```

The stored orders at `N=130,170` are reversed; the analyzer explicitly corrects them. The fresh results are:

| lineage | ratio | ratio SE | fixed residual z |
|---|---:|---:|---:|
| `65 -> 130` | -0.31382 | 0.0908 | +0.114 |
| `85 -> 170` | -0.34095 | 0.1118 | -0.150 |

The covariance-aware joint residual score is

\[
\chi^2=0.03445\quad(2\ \mathrm{dof}).
\]

This is strong independent evidence for the joint angular-phase and `13/8` magnitude relation at these two exact lineages.

### Harmonic boundary of the doubling result

For square-lattice harmonics,

\[
\cos[4m(\theta+\pi/4)]=(-1)^m\cos(4m\theta).
\]

The observed sign flip selects the odd-`m` class (`cos4`, `cos12`, `cos20`, ...) and rejects dominance by the even-`m` class (`cos8`, `cos16`, ...). It does not uniquely distinguish H4 from H12. H12 must therefore be retained as an explicit Stage-2 falsification alternative.

### What P37 does not establish

P37 is measured at one fixed thermal coordinate. It does not yet establish:

- the same ratio for the thermal-even matching-odd full-curve projector;
- the slope ratio `2^(3/8)`;
- the induced root-gap ratio `-1/4`;
- absence of finely aligned logarithmic or subleading odd-`m` sectors;
- a unique LCFT operator.

Issue #49 freezes these full-curve and root-ratio tests.

## 4. Threshold-rank infrastructure passes correctness and throughput gates

The C++17/OpenMP P33 engine implements the frozen bidirectional rank convention and stores integer `K_minus`/`K_plus` histograms by orientation and batch. It passes:

- all 120 permutations of the primitive Gaussian `N=5` torus;
- bin-for-bin Python-oracle comparisons on shared counters;
- byte-identical one-thread and two-thread output;
- `K_minus <= K_plus` on every sample.

The Huawei pilot processed 10,000,000 paired permutations per size for all five sizes in 65.2 seconds. The retained histograms reconstruct the full matching curves, analytic slopes, and roots.

The thermal-parity result is positive: at `p_ref +/- 0.001`, the central-reflection even orientation component persists, while the odd component agrees with `0.001*Delta M'`. The fixed-p signal is therefore not produced by a small common error in the thermal coordinate.

However, the radial root-gap gate does not pass. `N^2 Delta p` drifts, and a constant trained on `N=65,85,130` gives held-out chi-square `7.40` at `N=145,170`.

## 5. P35 closes local root conversion, not the radial exponent

From the same P33 histograms,

\[
C_N=-\frac{\Delta p_N^*\,\overline{M'_N}}{\Delta M_N}
\]

lies between `0.99984` and `1.00031`. Direct and linearized root gaps agree, and

\[
B_N=N^{-3/8}\overline{M'_N}
\]

changes only from approximately 1.7514 to 1.7462.

This excludes nonlinear root conversion or a large orientation-slope asymmetry as the cause of the radial drift. It is not an independent exponent measurement: `Delta M`, both roots, and `M'` come from the same finite curves. The amplitude relation

\[
A_p(N)=A_M(N)/B_N
\]

must still pass a cross-size, covariance-aware held-out test.

Issue #49 adds a sharper no-fit root prediction. If `M'_{2N}/M'_N=2^{3/8}`, then

\[
\boxed{
\frac{\Delta p^*_{2N}}{\Delta p^*_N}=-\frac14.
}
\]

## 6. Euler/motif controls are useful but not production-ready

The exact Euler identity and fixed-K hypergeometric centering pass on axis `L=3`, diamond `L=2`, and primitive Gaussian `N=5`. Pilot-frozen single-geometry variance reductions are:

```text
L=8   2.319x
L=12  1.852x
L=16  1.665x
```

The multiple-size `>=2x` gate remains closed. Issue #40 therefore targets the actual paired orientation difference using exact zero-mean differences of equal-multiplicity edge, face, and motif counts.

## 7. Reproducibility and covariance corrections

### 7.1 Fixed-p provenance

P31 metadata did not identify its executed five-design source by a clean commit. P37 now supplies a clean-source, independent 100-million replay at `N=65,85,130,170`. These four sizes should not be rerun solely for provenance wording.

Remaining fixed-p work under #39 is:

- executable SHA-256 and dirty-tree record for the P37 binary;
- a clean-source `N=145` replay;
- a per-seed compatibility table without premature pooling.

### 7.2 Threshold-rank provenance and implicit cross-size coupling

P33 records a working-tree source. Its `counter_permutation()` stream derives from `(seed, replica)` without mixing `N`; aligned counters therefore introduce implicit cross-size coupling.

Before the production replay, choose one policy:

1. domain-separate sizes by `N`/engine tag or disjoint seed domains; or
2. retain deliberate coupling and use the measured full cross-size covariance.

PR #46 adds an audit that reconstructs `Delta M`, slopes, root jackknife pseudo-values, `A_M`, `B`, and `A_p` from the existing aligned P33 batches. Its output must be committed before the production RNG choice. The output contract should also retain cross-orientation rank moments.

## 8. Frozen next experiments

### 8.1 Audit existing P33 covariance

Run PR #46 on the committed 10-million P33 histogram and compare full-covariance versus diagonal scores for both `A_M` and `A_p`. Do not replace the earlier diagonal result; report both.

### 8.2 Clean production threshold ranks

Complete #39/#26 with clean source and binary hashes, an explicit RNG-domain policy, at least 100,000,000 paired permutations per size, at least 100 batches, integer histograms, cross-orientation moments, and full jackknife covariance.

### 8.3 Full-curve doubling and root-ratio test

Issue #49 freezes thermal levels

```text
u = 0, 0.025, 0.050
```

and tests, in exact lineage order,

\[
X^{\rm even}_{4,2N}(u)+2^{-13/8}X^{\rm even}_{4,N}(u)=0,
\]

\[
\frac{\overline M'_{2N}}{\overline M'_N}=2^{3/8},
\qquad
\frac{\Delta p^*_{2N}}{\Delta p^*_N}=-\frac14.
\]

Use delete-one jackknife pseudo-values for roots and nonlinear ratios.

### 8.4 General radial root-gap challenge

Using clean rank data, compare on training sizes `N=65,85,130` and held-out `N=145,170`:

1. the H4-derived `N^-2` prediction using measured `B_N`;
2. one training-selected inverse-power correction;
3. a logarithmic correction;
4. a free root exponent trained only on the declared training set.

### 8.5 Stage-2 angular-radial challenge

Freeze:

```text
N=185  (13,4)/(11,8)   training extension
N=221  (14,5)/(11,10)  held out
N=265  (16,3)/(12,11)  held out
```

Keep the Stage-1 model order and add one explicit H4+H12 alternative before any new output is viewed. Retain H8 as an independent bound. The four-angle `N=1105` experiment remains gated on the two new held-out sizes and #49.

### 8.6 Paired local controls

Issue #40 tests exact zero-mean differences of NN-edge, face, and motif-orbit counts on the same-N pairs. Its gate is evaluated against the unadjusted paired `Delta M` estimator on fresh data.

### 8.7 Same-modulus `kappa3`

Square-bond and triangular-rhombus sequences have different torus moduli and remain method controls. The next universality test must match the macroscopic modulus and observable normalization. `-5/3` remains a candidate, not a result.

## 9. Operator interpretation boundary

The finite-size evidence now supports:

- a nonzero matching-odd orientation sector;
- independent-seed reproducibility;
- held-out preference for `Delta cos4 theta N^-13/8`;
- a no-fit `N -> 2N`, `theta -> theta+pi/4` magnitude and sign relation;
- odd harmonic parity under the `pi/4` rotation.

This makes a spin-4 field of total dimension

\[
x=21/4
\]

a concrete operator candidate. It is not yet a unique identification. A level-four thermal-family descendant or quasiprimary, possibly with a logarithmic partner in the `c=0` theory, remains conjectural until:

- #49 passes on the full-curve projector and root ratio;
- #36 bounds H12 and other odd harmonics;
- the selection-rule derivation excludes lower matching-odd sectors;
- modulus dependence is consistent.

## 10. Hardware decision

The fixed-p work is already fast on the 16-core CPU. The threshold-rank CPU engine is also efficient enough for the clean production replay. GPU work is justified only after clean CPU provenance, exact CPU/GPU equality, and an end-to-end information-per-wall-time benchmark.

Large-Pell scans and `N=1105` remain forbidden as first GPU workloads.

## 11. Decision

The project is now in mechanism-confirmation rather than effect-confirmation:

- **confirmed on the tested finite sizes:** nonzero same-N effect, predicted sign, independent seed, Stage-1 held-out H4 success, and the two-lineage no-fit doubling relation;
- **confirmed locally:** thermal-center robustness and linear residual-to-root conversion;
- **still open:** full-curve doubling, root ratio `-1/4`, unique asymptotic exponent, H4 versus H12, threshold-rank provenance/covariance, multiple-size variance gate, and the `x=21/4` operator identity.

Immediate execution order:

1. PR #46 existing-P33 covariance audit;
2. #39/#26 clean production ranks;
3. #49 full-curve doubling and root-ratio test;
4. #35 general radial root challenge;
5. #36 Stage 2 with H12;
6. #40 paired motif controls in parallel;
7. #25 same-modulus `kappa3`;
8. #12 operator selection only after the numerical gates pass.
