# Synthetic model red-team at the current Gaussian-orientation precision

This is a design-power calculation, not a fit to target data. It mirrors the frozen
P32 split (`65,85,130` train; `145,170` held out), uses the pooled P31 standard
errors, and preserves the exact Gaussian-pair H4/H12 angular columns.

Two selectors are reported: raw held-out chi-square (closest to the current fixed-model
score) and predictive deviance (chi-square plus forecast-volume and parameter penalties).
A model can remain statistically acceptable without being selected.

## Design

| N | split | SE | delta cos4 | delta cos12 |
|---:|:---|---:|---:|---:|
| 65 | train | 7.07e-05 | 1.363314 | -0.920117 |
| 85 | train | 7.83e-05 | 1.594464 | -0.264211 |
| 130 | train | 8.95e-05 | 1.363314 | -0.920117 |
| 145 | heldout | 7.71e-05 | 1.917717 | 1.303807 |
| 170 | heldout | 9.2e-05 | 1.594464 | -0.264211 |

## Noise scale 1x current SE; admixture 0.5

| hidden truth | correct: chi-square | correct: predictive deviance | H4 not rejected |
|:---|---:|---:|---:|
| H4 | 30.3% | 89.0% | 95.0% |
| H12 | 50.1% | 93.4% | 0.0% |
| H4+H12 | 92.9% | 92.4% | 0.0% |
| two_radial_powers | 37.1% | 14.1% | 79.5% |
| log_Jordan | 23.6% | 25.1% | 0.0% |
| ordinary_correction | 21.1% | 26.0% | 67.7% |

Predictive-deviance confusion matrix (row truth, column selected):

| truth \ selected | H4 | H12 | H4+H12 | two_radial_powers | log_Jordan | ordinary_correction |
|:---|---:|---:|---:|---:|---:|---:|
| H4 | 89.0% | 0.0% | 4.8% | 3.0% | 0.7% | 2.5% |
| H12 | 0.0% | 93.4% | 6.5% | 0.0% | 0.0% | 0.1% |
| H4+H12 | 0.0% | 0.0% | 92.4% | 1.7% | 1.4% | 4.6% |
| two_radial_powers | 55.8% | 0.0% | 11.7% | 14.1% | 4.0% | 14.4% |
| log_Jordan | 0.0% | 0.0% | 3.6% | 34.9% | 25.1% | 36.4% |
| ordinary_correction | 38.6% | 0.0% | 15.5% | 14.4% | 5.5% | 26.0% |

## Noise scale 0.5x current SE; admixture 0.5

| hidden truth | correct: chi-square | correct: predictive deviance | H4 not rejected |
|:---|---:|---:|---:|
| H4 | 30.3% | 89.6% | 94.9% |
| H12 | 50.1% | 93.1% | 0.0% |
| H4+H12 | 99.9% | 99.9% | 0.0% |
| two_radial_powers | 47.4% | 32.6% | 31.4% |
| log_Jordan | 46.6% | 47.7% | 0.0% |
| ordinary_correction | 37.8% | 47.2% | 10.7% |

Predictive-deviance confusion matrix (row truth, column selected):

| truth \ selected | H4 | H12 | H4+H12 | two_radial_powers | log_Jordan | ordinary_correction |
|:---|---:|---:|---:|---:|---:|---:|
| H4 | 89.6% | 0.0% | 4.7% | 2.5% | 0.7% | 2.3% |
| H12 | 0.0% | 93.1% | 6.9% | 0.0% | 0.0% | 0.0% |
| H4+H12 | 0.0% | 0.0% | 99.9% | 0.0% | 0.1% | 0.1% |
| two_radial_powers | 11.2% | 0.0% | 14.7% | 32.6% | 10.3% | 31.3% |
| log_Jordan | 0.0% | 0.0% | 0.0% | 30.6% | 47.7% | 21.8% |
| ordinary_correction | 4.3% | 0.0% | 19.3% | 17.9% | 11.3% | 47.2% |

## Interpretation guardrails

- Pure H12 is amplitude-normalized to the H4 signal norm; mixed alternatives add a
  component whose covariance-weighted norm is the declared admixture fraction of H4.
- `two_radial_powers` means H4 angular structure with N^-13/8 and N^-9/8 terms.
  `ordinary_correction` uses N^-13/8 and N^-21/8; `log_Jordan` uses
  N^-13/8 times 1 and log(N/Nref).
- Only two held-out points exist. Low correct-selection power or high H4 non-rejection
  is therefore a property of the current design, not evidence that the mechanisms are equal.
- Re-run after target covariance is available; do not tune admixture after seeing targets.
