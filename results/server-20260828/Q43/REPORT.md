# Q43 — Joint angular–radial harmonic challenge on C05 histograms

Issue 36. **Existing C05 data only.** Seed `0xC0100001`, **2e6** samples,
cross-channel Newman-Ziff `ΔM`. Primary exponent frozen at `13/8`.
Train N=65,85,130; hold out N=145,170. N=145 is not dropped.
Either/cross are not two replications. Second-stage N=185/221/265 MC skipped.

## Design columns (exact rationals)

| N | split | Δcos4 | Δcos8 | ΔM ± SE | N^{13/8} ΔM ± SE |
|---|---|---|---|---|---|
| 65 | train | 1152/845 | 3838464/3570125 | 0.000218787±0.000403792 | 0.193199±0.356567 |
| 85 | train | 2304/1445 | -10386432/10440125 | 0.00164377±0.000370388 | 2.24464±0.505781 |
| 130 | train | 1152/845 | -3838464/3570125 | 9.67594e-05±0.000460184 | 0.263543±1.2534 |
| 145 | holdout | 8064/4205 | -9257472/88410125 | 0.000226332±0.00052912 | 0.736156±1.72099 |
| 170 | holdout | 2304/1445 | 10386432/10440125 | 0.00091405±0.000456834 | 3.8499±1.92414 |

Selected power-correction `ω` (training only): **4.0**
(conditioning warning: True; ω=4 on N=65..130 is ill-conditioned.)
Free `α` (training only): **0.5**
(search interval [0.5, 3.0]; hit_lower_bound=True. A bound hit is not a physical exponent.)

## Model scores (covariance-aware χ² on ΔM)

| model | params | train χ² | holdout χ² | cond | notes |
|---|---|---|---|---|---|
| `fixed_13_8_cos4` | A4=0.632906±0.197105 | 9.72265 | 2.25669 | 1 |  |
| `fixed_13_8_cos4_power` | A4=2.24499±0.707965, B_power=-8458.72±3567.87 | 4.10194 | 2.78917 | 3.27661e+08 | ω=4.0 |
| `fixed_13_8_cos4_log_amplitude` | A4=-8.42423±4.78839, B_log=2.10151±1.1101 | 6.13889 | 3.38423 | 12191.1 |  |
| `fixed_13_8_cos4_plus_cos8` | A4=0.765552±0.202737, A8=-0.776856±0.277923 | 1.90939 | 3.34958 | 2.21577 | A8 bounded |
| `free_alpha_cos4` | A4=0.00510676±0.00147297 | 8.01329 | 1.63228 | 1 | α=0.5 |

**Raw lowest holdout χ²:** `free_alpha_cos4`
**Well-conditioned holdout winner (cond<1e4):** `fixed_13_8_cos4`
**Reported winner:** `fixed_13_8_cos4`
H8/log overturn H4 on holdout (Δχ²>1): **False**

H4 holdout residuals have the predicted sign at both N=145 and N=170.
Holdout χ² values are all O(1–3) on two noisy points; 2e6 does not distinguish
models sharply. `free_alpha` hits the preregistered lower bound α=0.5 and is
not treated as a physical winner. Ill-conditioned power/log companions are retained
as failed candidates.

### Held-out signed residuals (ΔM)

| model | N | true | pred | signed residual | z |
|---|---|---|---|---|---|
| `fixed_13_8_cos4` | 145 | 0.000226332 | 0.000373165 | -0.000146832 | -0.277503 |
| `fixed_13_8_cos4` | 170 | 0.00091405 | 0.000239593 | 0.000674457 | 1.47637 |
| `fixed_13_8_cos4_power` | 145 | 0.000226332 | 0.00108645 | -0.000860117 | -1.62556 |
| `fixed_13_8_cos4_power` | 170 | 0.00091405 | 0.000739063 | 0.000174986 | 0.383042 |
| `fixed_13_8_cos4_log_amplitude` | 145 | 0.000226332 | 0.00119951 | -0.000973176 | -1.83924 |
| `fixed_13_8_cos4_log_amplitude` | 170 | 0.00091405 | 0.000896697 | 1.7353e-05 | 0.0379855 |
| `fixed_13_8_cos4_plus_cos8` | 145 | 0.000226332 | 0.000476383 | -0.000250051 | -0.472579 |
| `fixed_13_8_cos4_plus_cos8` | 170 | 0.00091405 | 0.000106314 | 0.000807736 | 1.76812 |
| `free_alpha_cos4` | 145 | 0.000226332 | 0.000813291 | -0.000586959 | -1.10931 |
| `free_alpha_cos4` | 170 | 0.00091405 | 0.000624505 | 0.000289545 | 0.633809 |

## Is A8 resolved?

- status: **bounded**
- A8 = -0.776856 ± 0.277923 (z=-2.79522)
- corr(A4,A8) = -0.234071
- condition number = 2.21577
- identifiable design: True
- holdout χ² vs H4: 3.34958 vs 2.25669; improves=False
- rule: unidentifiable if cond>1e6 or |corr(A4,A8)|>0.98; resolved only if identifiable AND |A8|/se>=2 AND holdout chi2 beats H4; otherwise bounded (in-sample coefficient, not holdout-confirmed)

## Huawei published A4 (separate table, not pooled)

| N | Huawei A4 | C05 A_M (NZ M) | C05−Huawei | z (independent SEs) |
|---|---|---|---|---|
| 65 | 0.65±0.109 | 0.141713±0.261545 | -0.508287 | -1.79386 |
| 85 | 0.651±0.136 | 1.40777±0.317211 | 0.756769 | 2.19267 |
| 145 | 0.557±0.331 | 0.383871±0.897414 | -0.173129 | -0.181001 |

Huawei published A4 is a separate comparison table; never pooled with C05 as if same seed.

C05 `A_M` is the Newman-Ziff matching function; Huawei A4 is the fixed-p
wrapping difference. Different observable, different seed, different sample count.

## Provenance

- C05 source commit: `f89191a5468ec8417bbeac373334cecb6b5833a7`
- analysis branch: `agent/q42-q43-closure-harmonics`
- wall_time_s: 1.2042782306671143
- machine: cursor Linux 8 cores / 16 GB, no GPU

