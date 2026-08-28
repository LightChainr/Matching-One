# P48 retrospective projector score

**Status: retrospective planning only.** The P33 histogram source existed before P48 was frozen and its recorded source is a working tree. These numbers must not be presented as confirmatory or preregistered evidence.

## Decision

**Parity pattern supported but four-power joint not passed.** All four projected channels have the expected parity-resolved signal pattern and their frozen nonzero-amplitude predictions improve substantially on zero. However, the heldout `P4_S_prime` scaled amplitude drifts upward and fails its pure `N^-5/4` law; therefore the four-law package does not pass as a joint claim.

## Frozen train/heldout score

The amplitude of each law was fit using only `N=65,85,130`; `N=145,170` were scored as heldout. Prediction uncertainty includes amplitude-estimation uncertainty and the synchronized train/heldout cross covariance.

| channel | power | train amplitude (SE) | heldout chi2 / 2 | zero chi2 / 2 | heldout z (145, 170) |
|---|---:|---:|---:|---:|---:|
| P4_S | 1 | -0.0081257951 (0.00214) | 1.1083 | 9.3443 | -0.958, -0.427 |
| P4_D | 13/8 | 0.26491444 (0.0556) | 5.5290 | 18.6337 | +1.538, +1.812 |
| P4_S_prime | 5/4 | 1.9434248 (0.0766) | 10.1908 | 142.0115 | +1.963, +2.535 |
| P4_D_prime | 5/8 | -0.024523828 (0.00488) | 0.4073 | 10.7923 | -0.594, +0.299 |

The covariance-aware four-channel heldout omnibus statistic is **11.8933 / 8 df** (zero-model statistic 246.8687 / 8 df). The omnibus statistic alone is not unusually large; the declared package nevertheless fails its conjunction criterion because `P4_S_prime` is inconsistent with its frozen law (chi-square 10.1908 / 2 df). These are retrospective diagnostics, not confirmatory p-values.

## Method

At the intrinsic center where the mean of the two orientation matching functions vanishes, the analyzer reconstructs the thermal-even `S` and matching-odd `D` sectors and their analytic first derivatives. It projects the orientation contrast by the exact Gaussian-integer `Delta cos(4 theta)`. The 100 aligned batch ids are deleted synchronously across every N and channel; jackknife pseudovalues then supply the full covariance matrix used by GLS and predictive residual scoring.

The next legitimate use is design: freeze one correction/log alternative motivated by the `P4_S_prime` drift, then test it on fresh independent counters. Do not refit alternatives on the heldout values reported here.
