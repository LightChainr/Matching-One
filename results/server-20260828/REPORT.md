# Huawei ARM server analysis, 2026-08-28

This directory records the first combined run of the finite-size audit,
exact torus checks, and the bounded Pell-paired Monte Carlo discovery kernel.
The starting repository commit was
`fa145e88ccc72042c8a1a365a16b6127eed0be50`.

## Environment

- Huawei Cloud HCE 2 container, Linux 5.10, AArch64
- 16 HiSilicon CPU cores, 30 GiB usable RAM, no GPU
- Python 3.9.9 with mpmath 1.3.0, PyYAML 6.0.3, and SymPy 1.14.0
- GCC 10.3.1 with OpenMP for the Monte Carlo kernel

The complete machine and dependency record is in `environment.txt`.

## Finite-size blind audit

The server completed all 54 preregistered grid jobs:

- decimal precision: 60, 100, and 160 digits;
- minimum training width: 5 through 10;
- rolling holdout length: 2, 3, and 4;
- correction models: `(4)`, `(4,6)`, `(4,6,8)`, `(4,6,8,10)`, and
  `(4,6,8,10,12)`.

There were zero failed jobs. Selection used only rolling folds whose test data
ended at or before width 18; widths 19--21 were hidden until final scoring.
The selected configuration was:

```text
model powers:       4,6,8,10,12
minimum width:      8
rolling holdout:    2
precision:          160 digits
validation RMSE:    8.7900300416091e-12
validation span:    3.8416220260837e-11
```

The frozen fit then predicted widths 19--21 with:

```text
final-tail RMSE:          1.3775861250986e-11
maximum absolute error:   1.9374870639209e-11
training-only intercept:  0.59274605094603206266439366806726549
```

This is an out-of-sample finite-width prediction result, not a statistical
confidence interval for the infinite-lattice threshold. The full selection
record is `issue-5-summary.json`; all 54 raw JSON files, logs, and the run
manifest are under `issue-5-grid/`.

The final summary was regenerated after fixing a precision-ordering defect:
the selected 160-digit precision is now activated before decimal CSV strings
are parsed into `mp.mpf`. Model selection and the scientific conclusion are
unchanged. The corrected intercept is higher by `2.6933e-16`; repository
history preserves the superseded output.

## Exact matching checks

After replacing Python-version-specific population counting, the reference
enumerator runs on the server's Python 3.9.9. Exact weighted checks at
`p=0.37` reproduced the cluster-side and wrapping-side matching functions:

```text
axis L=3:     difference = -2.1084395886461e-81
diamond L=2:  difference = -3.3735033418338e-80
```

Exact polynomial enumeration also reproduced the recorded physical roots:

```text
axis L=4:     0.5906721123310282968959020114395128696211
diamond L=3:  0.5942523211685686997053812881560658246685
```

The machine-readable outputs are in `exact/`.

## Pell-paired fixed-p Monte Carlo

The OpenMP discovery kernel used shared counter-based random numbers for the
axis and diamond geometries at three probabilities centered on
`p_ref=0.59274605` with step `0.001`. Independent replicas require no
autocorrelation correction.

| Pell pair | replicas | axis linear root | diamond linear root | diamond-axis gap |
|---|---:|---:|---:|---:|
| `(17,12)` | 1,000,000 | `0.59277769 +/- 5.39e-5` | `0.59271409 +/- 6.25e-5` | `-6.36e-5 +/- 5.98e-5` |
| `(41,29)` | 500,000 | `0.59273997 +/- 3.73e-5` | `0.59275555 +/- 4.08e-5` | `+1.56e-5 +/- 5.26e-5` |
| `(99,70)` | 100,000 | `0.59273994 +/- 4.49e-5` | `0.59277265 +/- 4.29e-5` | `+3.27e-5 +/- 6.06e-5` |

The gap changes sign and is below roughly one standard error for every pair.
Therefore this run does **not** resolve the proposed orientation-odd asymptotic
amplitude. It is a useful negative planning result: a three-point wrapping
indicator scan has insufficient variance reduction for the tiny `L^-4`
orientation gap, even though it is fast. A production follow-up should use the
full Newman--Ziff microcanonical curve and the newly implemented multi-channel
homology/control-variate machinery rather than increasing this fixed-p sample
count blindly.

The common-random-number coupling reduced the root-gap standard error most for
`(17,12)` (about 28% relative to combining the two marginal errors as
independent), but the benefit was only a few percent for the two larger pairs.
The present kernel couples axis and diamond sites by deterministic vertex
index, not by an optimized spatial correspondence; a better geometric coupling
is another possible variance-reduction target.

The batch aggregates, provenance metadata, jackknife analyses, and flat metric
CSVs are under `pell/`.

### Small-Pell power-calibration follow-up

After the first campaign established that the large pairs were underpowered,
the omitted `(7,5)` pair was run twice with independent seeds and 5,000,000
replicas per run. The orientation root gaps were

```text
seed 202608280705:  +4.6267e-4 +/- 5.06e-5
seed 202608290705:  +4.0417e-4 +/- 5.10e-5
inverse-variance pool: +4.3365e-4 +/- 3.59e-5
```

The pooled signal is about 12 standard errors from zero and is compatible
between seeds. Repeating the first seed with the finite-difference step reduced
from `0.001` to `0.0005` changed the gap by only `7.9e-7`, far below its sampling
error.

Combining this discovery point with the exact `(3,2)` root gap gives a
two-point effective exponent of about `4.23`; the `(7,5)` gap times the fourth
power of the mean physical length is about `1.06`. Tiny systems are not an
asymptotic fit, but both numbers agree closely with the preregistered `L^-4`
power budget. This supports using `(7,5)` as a positive power-calibration
cross-check while retaining same-N tomography as the cleaner primary test.

## Validation and scope

- 15 Python/C++ integration and exact-regression tests passed on the server.
- Exact homology tests exhaust all 512 axis `L=3` configurations and all 256
  diamond `L=2` configurations.
- The new covariance estimator freezes weights on a pilot sample and evaluates
  on an independent sample; it is implemented and tested but was not applied
  to the fixed-p Pell files because axis and diamond finite-size observables do
  not have exactly equal means.
- The C++ kernel is deliberately a bounded discovery engine, not a completed
  Newman--Ziff implementation.

## Server-priorities v2 follow-up

The follow-up programs were run from base commit `fad82a0` with the new S0--S2
source files present in the working tree; those files are committed in this
pull request immediately after the run. The raw same-N metadata preserves the
base SHA and records this source-state qualification explicitly.

### S0: arbitrary integer-period topology

The general reference engine now represents a torus by any nonsingular integer
`2 x 2` period matrix. Closed lifted displacements are converted to period-basis
windings with exact adjugate/determinant arithmetic. Axis, diamond, Gaussian
circulant, and unimodular basis-change regressions pass; exhaustive det-5 tests
leave homology rank and either/cross wrapping invariant under the basis change.

### S1: same-N Gaussian orientation tomography

A C++17/OpenMP engine measured the prescribed same-N pairs with shared
counter-based random fields. A 2,000,000-replica pilot was followed by a frozen
30,000,000-replica production seed. For the `either` wrapping convention, the
production matching-function differences (first orientation minus second) are:

| N | representations | difference | batch SE | z | `N^(13/8) difference / delta cos(4 theta)` |
|---:|---|---:|---:|---:|---:|
| 65 | `(8,1)` / `(7,4)` | `+1.00377e-3` | `1.68e-4` | `5.96` | `0.650` |
| 85 | `(9,2)` / `(7,6)` | `+7.60333e-4` | `1.58e-4` | `4.80` | `0.651` |
| 145 | `(12,1)` / `(9,8)` | `+3.28667e-4` | `1.95e-4` | `1.68` | `0.557` |

All three signs agree with `delta cos(4 theta)`, and the first two normalized
amplitudes are nearly identical. The third is statistically weaker but
compatible. This is substantially cleaner than the large-Pell fixed-p scan and
supports continued same-N angular testing.

The production matching-even differences were much less stable (`z=4.19,
0.18,1.66`) and did not show the proposed clean dominance over the matching-odd
sector. Thus the working claim that the even spin-4 signal must be larger is
not supported by this run; it should not be used as a gate for interpreting the
resolved matching-function differences without further controls.

### S2: square-bond exact-threshold kappa3 control

The `p=1/2` score estimator was checked by exact enumeration and then run with
1,000,000 independent samples per size:

| L | method | kappa3 | jackknife SE |
|---:|---|---:|---:|
| 2 | exact | `-1.1705532693` | -- |
| 3 | exact | `-1.4555871991` | -- |
| 4 | Monte Carlo | `-1.56230` | `0.01483` |
| 6 | Monte Carlo | `-1.60727` | `0.02169` |
| 8 | Monte Carlo | `-1.57106` | `0.02584` |
| 12 | Monte Carlo | `-1.65012` | `0.03406` |
| 16 | Monte Carlo | `-1.68356` | `0.03909` |

A weighted fixed-`L^-3/2` fit on `L=4,6,8,12` predicts the held-out `L=16`
point within `1.55` standard errors. Fitting all five Monte Carlo sizes gives
an intercept near `-1.649 +/- 0.023` from statistical errors alone, so `-5/3`
is retained but neither established nor sharply tested. Union-jack and
other controls remained unimplemented in this first sequence; no universality
or rational-value claim is made from the square-bond sequence alone.

## Master-queue C01--C07 follow-up

- **C01/C02:** the five prescribed wrapping channels are now retained with
  full covariance. A fresh non-overlapping 2-million evaluation at `N=65,85`
  is compatible with the earlier 30-million matching-function signal, but the
  proposed matching-even dominance is not supported and held-out multi-angle
  confirmation remains necessary.
- **C03:** all five matching-odd channels are identical configuration by
  configuration. Frozen GLS weights give exactly `1.0x` variance reduction at
  `N=65,85`; this wrapping-only GPU gate fails for a structural reason, not
  lack of samples. That negative result is retained as provenance.
- **P34:** Euler and local-motif controls replace wrapping GLS. The identity
  `C_black-C_white=q+V-E+F0` holds on exhaustive tiny tori. Duplicate wrapping
  channels are rejected rather than GLS-combined. Pilot-frozen OLS on
  occupancy, NN edges, faces and short motifs yields about `2.34x` / `2.16x`
  variance reduction versus `q` at `N=65` / `N=85` on 800,000 fresh replicas.
  The redesigned `>=2x` GPU gate passes; GPU was not started. See `P34/`.
- **C04:** a triangular-site self-matching `p=1/2` control now passes exact
  `L=2,3,4` regressions. A 300,000-replica sequence through `L=32` validates
  the derivative pipeline but does not determine the correction exponent or
  establish `kappa3=-5/3`.
- **C05:** the threshold-rank reference freezes `K_minus/K_plus` conventions,
  reconstructs `M`, `M'`, and the root from integer histograms, and completed
  an axis `L=8`, 100,000-permutation ARM pilot. Its 83.85-second runtime makes
  a C++/GPU port preferable before multi-size production.
- **C07:** leakage-safe Stage A selected polynomial `F` degree 4 with
  `n_min=9`. Padé `[2/2]` was worse; both retained the positive, increasing
  signed errors at widths 19--21. Rational corrections did not cure the drift.

The final local and Python 3.9/AArch64 server regression suites each contain
42 tests; all 42 passed in both environments.

## Third-wave confirmation follow-up

- **P31:** a new independent 100-million-replica seed at each of
  `N=65,85,130,145,170` reproduced the prescribed sign at every size.  The new
  N=65 and N=85 signals are 16.0 and 11.2 standard errors.  Pooling seeds
  within size gives a common scaled amplitude `A4=0.7885 +/- 0.0352`, with
  chi-square `1.53/4` across all five sizes.
- **P32:** the frozen training/held-out challenge selected the simple fixed
  `13/8` model on held-out performance: chi-square `1.058` for the three
  held-out seed rows, versus `37.32` for a zero-effect model.  Power, logarithmic,
  H4+H8, and free-exponent alternatives did not improve prediction.  The H8
  coefficient is `-0.0345 +/- 0.0542` and is not resolved.
- **P33:** the C++ threshold-rank engine now retains orientation/batch
  `K_minus/K_plus` histograms and joint moments for all five pairs.  A
  10-million-per-size ARM pilot completed in 65.2 seconds.  Central amplitudes
  and root gaps have the predicted signs, and symmetric thermal offsets leave
  the central orientation amplitude intact.  Its historical constant raw
  `N^2 DeltaRoot` score was model-mismatched because angular leverage varies;
  the correctly normalized prospective score is reported under P45.
- **P34:** the preliminary axis run gave `2.319x`, `1.852x`, and `1.665x` at
  `L=8,12,16`.  The later same-N run added exact Gaussian controls and a richer
  microcanonical motif basis; on held-out replicas it gives about `2.34x` at
  `N=65` and `2.16x` at `N=85` versus the best single estimator.  This passes
  the multiple-size gate for individual-geometry matching estimates.  The
  orientation-difference OLS overfits its pilot and remains excluded.
- **P35:** same-batch amplitude closure gives
  `C=-Delta p* mean(M')/Delta M` between `0.99984` and `1.00031` at all five
  sizes, with delete-one-batch errors at the `2e-5--1.4e-4` level.  Direct and
  linearized root gaps agree; the remaining root-scaling failure is radial
  drift, not nonlinear root conversion.
- **P37:** the preregistered Gaussian doubling test predicts the parameter-free
  ratio `Delta M(2N)/Delta M(N)=-2^(-13/8)`.  A fresh 100-million-replica seed
  gives `-0.3138 +/- 0.0908` for `65->130` and `-0.3410 +/- 0.1118` for
  `85->170`.  The fixed-prediction residuals are only `+0.114` and `-0.150`
  standard errors, with joint chi-square `0.03445/2`.  This independently
  supports the combined spin-4 sign and exponent relation without fitting an
  amplitude.
- **P45:** a clean-source 100-million-per-size threshold-rank replay at
  `N=65,85` gives angular-normalized root amplitudes `0.42034 +/- 0.02157` and
  `0.39495 +/- 0.03078`.  With full cross-size jackknife covariance and the
  frozen source uncertainty, the preregistered `A_p=0.45101` prediction scores
  `chi-square=2.4267/2` and passes.  A free common amplitude is
  `0.41301 +/- 0.01924`; zero effect scores `461.3/2`.  Closure remains within
  `9e-5` of one, so the result tests the angular amplitude rather than root
  nonlinearity.
- **P48 retrospective:** synchronized cross-size jackknife reconstruction of
  the P33 histograms supports the alternating matching-parity channel pattern,
  but the four frozen powers do not pass as a conjunction.  In particular,
  held-out `P4[D'] ~ N^-5/8` scores `0.407/2`, whereas
  `P4[S'] ~ N^-5/4` scores `10.191/2` and shows upward drift.  This is planning
  evidence only because P33 predates the P48 protocol.
- **P50-A:** the prospective third Gaussian doubling lineage passes at
  `N=290`.  A fresh 500-million-replica run gives
  `Delta M=-1.60648e-4 +/- 4.0542e-5` versus the frozen
  `-1.37656e-4 +/- 2.4997e-5` target.  The combined residual is only
  `-0.483` standard errors (`chi-square=0.233/1`), while zero is excluded at
  `3.96` sampling standard errors.  This extends the no-fit `-2^-13/8`
  relation to a third exact lineage.

The expanded local suite contains 67 tests; all passed.  Detailed reports and
machine-readable outputs are under `P31/`, `P32-radial-challenge/`, `P33/`,
`P34/`, `P35-amplitude-closure/`, `P37-doubling/`, `P45-root-amplitude/`,
`P48-retrospective/`, `P50-third-lineage/`, and `C03-euler/`.
