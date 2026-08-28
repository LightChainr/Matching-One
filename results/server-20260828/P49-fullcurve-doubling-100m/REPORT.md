# Clean 100M full-curve norm-2 score for issues #48/#49

**Classification:** fresh-seed replication at geometries already present in the retrospective P33 source. No target parameter was fit.

## Decision

The full-curve thermal-even and root tests remain compatible but are driven to about 2.1 sigma by the 85->170 lineage. The raw asymptotic slope multiplier is sharply rejected at this precision, exposing a small finite-size correction. The four-channel P48 pure-power conjunction fails specifically and reproducibly in `P4[S']`; both preregistered correction models survive this same-geometry fresh-seed replication.

## P49 no-fit Gaussian doubling

| statistic | chi-square / df |
|---|---:|
| `P49_X_even_u=0.0` | 4.44798 / 2 |
| `P49_X_even_u=0.025` | 4.45985 / 2 |
| `P49_X_even_u=0.05` | 4.49571 / 2 |
| `P49_root_raw` | 4.48059 / 2 |
| `P49_root_finite_slope` | 4.44792 / 2 |
| `P49_slope` | 6412.89 / 2 |

Point ratios:

- `65->130`: slope=1.29398351 (target 1.29683955); root=-0.268717842 (raw target -0.25).
- `85->170`: slope=1.29437757 (target 1.29683955); root=-0.402402308 (raw target -0.25).

All `u={0,.025,.05}` coordinates were solved inside every delete-one replicate. The reported inverse-slope condition proxies and coordinates are retained in `analysis/score.json`.

## P48 derivative spectrum

The mathematically consistent child/parent factor for a projector normalized by each size's own `Delta cos(4 theta)` is **positive** `2^-alpha`. The negative H4 rotation factor belongs to the unnormalized exact-lineage contrast. Applying a negative factor to normalized `P4` double-counts the angular sign.

| normalized channel | positive-ratio chi-square / 2 | negative-artifact chi-square / 2 |
|---|---:|---:|
| `P4_S` | 0.0587199 | 143.649 |
| `P4_D` | 4.44798 | 175.282 |
| `P4_S_prime` | 39.9817 | 3363.73 |
| `P4_D_prime` | 4.68483 | 99.9921 |
| **joint** | **76.9665 / 8** | **5595.34 / 8** |

The normalized pure laws pass cleanly for `P4[S]`, are mildly strained for `P4[D]` and `P4[D']`, and fail for `P4[S']`.

## Frozen P48 S-prime models on fresh N=130/170 counters

| model (frozen order) | chi-square / 2 | marginal z at N=130,170 |
|---|---:|---:|
| `pure_N^-5/4` | 37.8868 | +3.872, +5.747 |
| `zero_effect` | 1044.57 | +24.324, +21.394 |
| `q2_even_scalar` | 1.78993 | -0.421, +0.653 |
| `rank2_jordan_log` | 0.676706 | -0.524, +0.083 |

The q=2 correction is the first surviving model in the frozen order. The Jordan-log adversary also survives and has a smaller descriptive chi-square, but these same geometries do not constitute a new angular holdout and should not be used to reverse the preregistered order.

## Provenance

- Parent N=65/85 and child N=130/170 use source commit `6d2d68a`, seed `2026104501`, counters `[5000000000,5100000000)`, 100 aligned batches, and 100,000,000 samples per orientation pair.
- Child runs were executed on Huawei DevEnv `f415a4bcbd9a438b85f5f29e4a507ea4` (AArch64, 16 vCPU); both stderr files are empty.
- Raw child histograms, moments, metadata, stdout and stderr are preserved under `raw/`.
