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
final-tail RMSE:          1.3775763431649e-11
maximum absolute error:   1.9374735274555e-11
training-only intercept:  0.59274605094603179333564156705784225
```

This is an out-of-sample finite-width prediction result, not a statistical
confidence interval for the infinite-lattice threshold. The full selection
record is `issue-5-summary.json`; all 54 raw JSON files, logs, and the run
manifest are under `issue-5-grid/`.

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
